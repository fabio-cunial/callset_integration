version 1.0


# Downloads BAM files up to a given coverage.
#
workflow HGSVC2Download {
    input {
        String sample_id
        Array[String] bam_addresses
        Int target_coverage
        String remote_dir
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
        remote_dir: "Root directory in the remote bucket. Every sample is stored in a subdirectory."
    }
    
    call HGSVC2DownloadImpl {
        input:
            sample_id = sample_id,
            bam_addresses = bam_addresses,
            target_coverage = target_coverage,
            remote_dir = remote_dir,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb
    }
    
    output {
    }
}


task HGSVC2DownloadImpl {
    input {
        String sample_id
        Array[String] bam_addresses
        Int target_coverage
        String remote_dir
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = (3*target_coverage)*2 + 512
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
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        BILLING_PROJECT="broad-firecloud-dsde-methods"
        
        LIST_FILE=~{write_lines(bam_addresses)}
        cat ${LIST_FILE}
        TARGET_N_BYTES=$(( ~{target_coverage} * 3000000000 * 2 ))
        touch tmp1.fastq
        while read ADDRESS; do
            if [[ ${ADDRESS} == gs://* ]]; then
                ${TIME_COMMAND} gsutil -u ${BILLING_PROJECT} cp ${ADDRESS} .
            else
                ${TIME_COMMAND} wget ${ADDRESS}
            fi
            FILE_NAME=$(basename ${ADDRESS})
            if [[ ${FILE_NAME} == *.bam ]]; then
                ${TIME_COMMAND} samtools fastq -@ ${N_THREADS} -n ${FILE_NAME} >> tmp1.fastq 
            elif [[ ${FILE_NAME} == *.fastq.gz ]]; then
                gunzip -c ${FILE_NAME} >> tmp1.fastq
            elif [[ ${FILE_NAME} == *.fastq ]]; then
                cat ${FILE_NAME} >> tmp1.fastq
            fi
            rm -f ${FILE_NAME}
            N_BYTES=$(wc -c < tmp1.fastq)
            if [ ${N_BYTES} -gt ${TARGET_N_BYTES} ]; then
                break
            fi
        done < ${LIST_FILE}
        head -c ${TARGET_N_BYTES} tmp1.fastq > tmp2.fastq
        rm -f tmp1.fastq
        N_ROWS=$(wc -l < tmp2.fastq)
        N_ROWS=$(( (${N_ROWS}/4)*4 ))
        head -n ${N_ROWS} tmp2.fastq | gzip > ~{sample_id}.fastq.gz
        rm -f tmp2.fastq
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.fastq.gz ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
