version 1.0


#
workflow TruvariIntersampleNaive {
    input {
        Array[File] intrasample_merged_vcf
        Array[File] intrasample_merged_vcf_tbi
    }
    parameter_meta {
    }

    call TruvariIntersampleNaiveImpl {
        input:
            intrasample_merged_vcf = intrasample_merged_vcf,
            intrasample_merged_vcf_tbi = intrasample_merged_vcf_tbi
    }
    
    output {
    	File bcftools_merged = TruvariIntersampleNaiveImpl.bcftools_merged
    	File bcftools_merged_idx = TruvariIntersampleNaiveImpl.bcftools_merged_idx
        File truvari_collapsed = TruvariIntersampleNaiveImpl.truvari_collapsed
    	File truvari_collapsed_idx = TruvariIntersampleNaiveImpl.truvari_collapsed_idx
    }
}


task TruvariIntersampleNaiveImpl {
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
        
        INPUT_FILES=~{sep=',' intrasample_merged_vcf}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for FILE in ${INPUT_FILES}; do
            ID=$(basename ${FILE} .vcf.gz)
            bcftools norm --multiallelics - --output-type z ${FILE} > ${ID}_normed.vcf.gz
            tabix ${ID}_normed.vcf.gz
            echo ${ID}_normed.vcf.gz >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > bcftools_merged.vcf.gz 
        tabix bcftools_merged.vcf.gz
        ${TIME_COMMAND} bcftools norm --multiallelics - --output-type z bcftools_merged.vcf.gz > bcftools_merged_normed.vcf.gz 
        tabix bcftools_merged_normed.vcf.gz
        ${TIME_COMMAND} truvari collapse -i bcftools_merged_normed.vcf.gz -c removed.vcf.gz \
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
        memory: "128GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
