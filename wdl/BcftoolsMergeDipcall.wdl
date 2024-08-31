version 1.0


#
workflow BcftoolsMergeDipcall {
    input {
        Array[String] sample_id
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        Int min_sv_length
    }
    parameter_meta {
    }
    
    call InterSampleMerge {
        input:
            sample_id = sample_id,
            input_vcf_gz = sample_vcf_gz,
            input_tbi = sample_tbi,
            min_sv_length = min_sv_length
    }
    
    output {
        File output_vcf_gz = InterSampleMerge.output_vcf_gz
        File output_tbi = InterSampleMerge.output_tbi
    }
}


task InterSampleMerge {
    input {
        Array[String] sample_id
        Array[File] input_vcf_gz
        Array[File] input_tbi
        Int min_sv_length
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int n_files = length(input_vcf_gz)
    
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
        SAMPLE_IDS=~{sep=',' sample_id}
        rm -f list.txt
        for i in $(seq 1 ~{n_files}); do
            # Enforcing the right sample name
            INPUT_FILE=$(echo ${INPUT_FILES} | cut -d , -f ${i})
            SAMPLE_ID=$(echo ${SAMPLE_IDS} | cut -d , -f ${i})
            echo ${SAMPLE_ID} > samples.txt
            bcftools reheader --samples samples.txt ${INPUT_FILE} > tmp1.vcf.gz
            tabix -f tmp1.vcf.gz
            # Removing multiallelic records
            bcftools norm --multiallelics - --output-type v tmp1.vcf.gz > tmp2.vcf
            rm -f tmp1.vcf.gz*
            # Keeping only long events
            java -cp ~{docker_dir} Dipcall2VCF tmp2.vcf ~{min_sv_length} ${SAMPLE_ID}.fixed.vcf
            rm -f tmp2.vcf
            bgzip ${SAMPLE_ID}.fixed.vcf
            tabix -f ${SAMPLE_ID}.fixed.vcf.gz
            echo ${SAMPLE_ID}.fixed.vcf.gz >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        # Removing multiallelic records
        bcftools norm --multiallelics - --output-type v tmp1.vcf.gz > tmp2.vcf
        rm -f tmp1.vcf.gz*
        # Keeping only long events
        java -cp ~{docker_dir} Dipcall2VCF tmp2.vcf ~{min_sv_length} merged.vcf
        rm -f tmp2.vcf
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
        cpu: 8
        memory: "32GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
