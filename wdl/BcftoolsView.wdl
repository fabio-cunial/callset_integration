version 1.0


# Tries to implement $bcftools view$ as fast as possible on a large VCF.
#
workflow BcftoolsView {
    input {
        File vcf_gz
        File vcf_tbi
        File regions_bed
    }
    parameter_meta {
        vcf_gz: "Assumed to be sorted"
    }
    
    call BcftoolsViewImpl {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            regions_bed = regions_bed
    }
    
    output {
        File output_vcf_gz = BcftoolsViewImpl.output_vcf_gz
        File output_tbi = BcftoolsViewImpl.output_tbi
    }
}


# Performance on the 1074 AoU 8x VCF:
# COMMAND           RUNTIME     N_CPUS      MAX_RSS
# bcftools view     
#
task BcftoolsViewImpl {
    input {
        File vcf_gz
        File vcf_tbi
        File regions_bed
    }
    parameter_meta {
        vcf_gz: "Assumed to be sorted"
    }
    
    Int ram_size_gb = 2*ceil(size(vcf_gz, "GB"))
    Int disk_size_gb = 4*ceil(size(vcf_gz, "GB"))
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --compression-level 1 --output-type z --regions-file ~{regions_bed} ~{vcf_gz} > out.vcf.gz
        tabix -f out.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " SSD"
        preemptible: 0
    }
}
