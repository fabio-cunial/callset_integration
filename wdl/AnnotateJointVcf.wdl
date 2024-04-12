version 1.0


# 
#
workflow AnnotateJointVcf {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
    }
    parameter_meta {
    }

    call AnnotateJointVcfImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi
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
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int disk_size_gb = 10*( ceil(size(input_vcf_gz,"GB")) )

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        export BCFTOOLS_PLUGINS="~{docker_dir}/bcftools-1.19/plugins"
        
        
        # Reming fields from INFO
        #
        #
        
        # Adding SVTYPE, SVLEN.
        ${TIME_COMMAND} truvari anno svinfo --minsize 0 --output tmp1.vcf.gz ~{input_vcf_gz}
        rm -f ~{input_vcf_gz}
        
        # Adding GT count: UNK, REF, HET, HOM.
        ${TIME_COMMAND} truvari anno gtcnt --output tmp2.vcf.gz tmp1.vcf.gz
        rm -f tmp1.vcf.gz
        
        # Adding all bcftools tags
        ${TIME_COMMAND} bcftools +fill-tags tmp2.vcf.gz -Oz -o annotated.vcf.gz -- -t all
        rm -f tmp2.vcf.gz
        tabix -f annotated.vcf.gz
    >>>

    output {
        File annotated_vcf_gz = work_dir + "/annotated.vcf.gz"
        File annotated_tbi = work_dir + "/annotated.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}