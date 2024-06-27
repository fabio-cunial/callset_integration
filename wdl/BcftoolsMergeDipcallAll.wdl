version 1.0


# Merges the entire raw dipcall files, which include SVs and SNPs.
#
workflow BcftoolsMergeDipcallAll {
    input {
        Array[File] sample_vcf_gz
    }
    parameter_meta {
    }
    
    call BcftoolsMergeDipcallAllImpl {
        input:
            input_vcf_gz = sample_vcf_gz
    }
    
    output {
        File output_vcf_gz = BcftoolsMergeDipcallAllImpl.output_vcf_gz
        File output_tbi = BcftoolsMergeDipcallAllImpl.output_tbi
    }
}


task BcftoolsMergeDipcallAllImpl {
    input {
        Array[File] input_vcf_gz
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        SAMPLE_ID="0"
        for INPUT_FILE in ${INPUT_FILES}; do
            SAMPLE_ID=$(( ${SAMPLE_ID} + 1 ))
            tabix -f ${INPUT_FILE}
            bcftools norm --multiallelics - --output-type v ${INPUT_FILE} > ${SAMPLE_ID}.vcf
            bgzip ${SAMPLE_ID}.vcf
            tabix -f ${SAMPLE_ID}.vcf.gz
            echo ${SAMPLE_ID}.vcf.gz >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # Removing multiallelic records one last time
        bcftools norm --multiallelics - --output-type v tmp1.vcf.gz > merged.vcf
        rm -f tmp1.vcf.gz*
        bgzip merged.vcf
        tabix -f merged.vcf.gz
        ls -laht; tree
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "64GB"
        disks: "local-disk 500 HDD"
        preemptible: 0
    }
}
