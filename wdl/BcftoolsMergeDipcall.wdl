version 1.0


#
workflow BcftoolsMergeDipcall {
    input {
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
    }
    parameter_meta {
    }
    
    call InterSampleMerge {
        input:
            input_vcf_gz = sample_vcf_gz,
            input_tbi = sample_tbi
    }
    
    output {
        File output_vcf_gz = InterSampleMerge.output_vcf_gz
        File output_tbi = InterSampleMerge.output_tbi
    }
}


task InterSampleMerge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
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
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > merge.vcf.gz
        tabix merge.vcf.gz
        ls -laht; tree
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merge.vcf.gz"
        File output_tbi = work_dir + "/merge.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: "32GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
