version 1.0


# Ensures that the final inter-sample VCF has a simple but consistent format.
#
workflow AnnotateJointVcf {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File new_header
    }
    parameter_meta {
        new_header: "Custom header to force on the file. Assumed not to include the annotations that are added by this task. Must not include the '#CHROM' line."
    }

    call AnnotateJointVcfImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            new_header = new_header
    }
    
    output {
        File annotated_vcf_gz = AnnotateJointVcfImpl.annotated_vcf_gz
        File annotated_tbi = AnnotateJointVcfImpl.annotated_tbi
    }
}


task AnnotateJointVcfImpl {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File new_header
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int disk_size_gb = 200*( ceil(size(input_vcf_gz,"GB")) )

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.19/plugins"
        
        gunzip -c ~{input_vcf_gz} > original.vcf

        ${TIME_COMMAND} java -cp ~{docker_dir} CleanIntersampleVcf original.vcf body.txt CollapseId NumCollapsed NumConsolidated
        cat ~{new_header} body.txt > tmp1.vcf
        rm -f body.txt
        head -n 500 tmp1.vcf
        bgzip tmp1.vcf
        tabix -f tmp1.vcf.gz
        
        # Adding SVTYPE, SVLEN.
        ${TIME_COMMAND} truvari anno svinfo --minsize 0 --output tmp2.vcf.gz tmp1.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # Adding GT count: UNK, REF, HET, HOM.
        ${TIME_COMMAND} truvari anno gtcnt --output tmp3.vcf.gz tmp2.vcf.gz
        tabix -f tmp3.vcf.gz
        rm -f tmp2.vcf.gz*
        
        # Adding all bcftools tags
        ${TIME_COMMAND} bcftools +fill-tags tmp3.vcf.gz -Oz -o annotated.vcf.gz -- -t all
        tabix -f annotated.vcf.gz
        rm -f tmp3.vcf.gz*
    >>>

    output {
        File annotated_vcf_gz = work_dir + "/annotated.vcf.gz"
        File annotated_tbi = work_dir + "/annotated.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}