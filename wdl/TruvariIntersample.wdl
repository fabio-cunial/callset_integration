version 1.0


#
workflow TruvariIntersample {
    input {
        Array[File] intrasample_merged_vcf
        Array[File] intrasample_merged_vcf_tbi
    }
    parameter_meta {
    }

    call TruvariIntersampleImpl {
        input:
            intrasample_merged_vcf = intrasample_merged_vcf,
            intrasample_merged_vcf_tbi = intrasample_merged_vcf_tbi
    }
    
    output {
    	File bcftools_merged = TruvariIntersampleImpl.bcftools_merged
    	File bcftools_merged_idx = TruvariIntersampleImpl.bcftools_merged_idx
        File truvari_collapsed = TruvariIntersampleImpl.truvari_collapsed
    	File truvari_collapsed_idx = TruvariIntersampleImpl.truvari_collapsed_idx
    }
}


task TruvariIntersampleImpl {
    input {
        Array[File] intrasample_merged_vcf
        Array[File] intrasample_merged_vcf_tbi
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(intrasample_merged_vcf,"GB"))) + 50
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        INPUT_FILES=~{sep='-' intrasample_merged_vcf}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr '-' ' ')
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none ${INPUT_FILES} --output-type z > bcftools_merged.vcf.gz 
        tabix bcftools_merged.vcf.gz
        ${TIME_COMMAND} truvari collapse -i bcftools_merged.vcf.gz -c removed.vcf.gz \
            --sizemin 0 --sizemax 1000000 -k common --gt all \
            | bcftools sort -m 2G --output-type z > truvari_collapsed.vcf.gz
        tabix truvari_collapsed.vcf.gz
    >>>
    
    output {
        File bcftools_merged = "~{work_dir}/bcftools_merged.vcf.gz"
        File bcftools_merged_idx = "~{work_dir}/bcftools_merged.vcf.gz.tbi"
        File truvari_collapsed = "~{work_dir}/truvari_collapsed.vcf.gz"
        File truvari_collapsed_idx = "~{work_dir}/truvari_collapsed.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "32GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
