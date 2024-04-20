version 1.0


# Remark: the workflow normalizes multiallelic sites and then filters SVs.
#
workflow Dipcall2SVs {
    input {
        File input_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
    }
    
    call Dipcall2SVsImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            min_sv_length = min_sv_length
    }
    
    output {
        File sv_vcf_gz = Dipcall2SVsImpl.sv_vcf_gz
        File sv_tbi = Dipcall2SVsImpl.sv_tbi
    }
}


task Dipcall2SVsImpl {
    input {
        File input_vcf_gz
        Int min_sv_length
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(input_vcf_gz,"GB"))
    Int ram_size_gb = 4
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
        
        bcftools norm --multiallelics - --output-type v ~{input_vcf_gz} > input.vcf
        rm -f ~{input_vcf_gz}
        ${TIME_COMMAND} java -cp ~{docker_dir} Dipcall2VCF input.vcf ~{min_sv_length} sv.vcf
        ${TIME_COMMAND} bgzip sv.vcf
        tabix -f sv.vcf.gz
    >>>
    
    output {
        File sv_vcf_gz = work_dir + "/sv.vcf.gz"
        File sv_tbi = work_dir + "/sv.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
