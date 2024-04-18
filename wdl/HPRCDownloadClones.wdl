version 1.0


# Downloads BAM files up to a given coverage, randomizing both the order in
# which the BAMs are downloaded, and the order of the reads in the
# concatenation of all downloaded BAMs.
#
workflow HPRCDownloadClones {
    input {
        String sample_id_clone
        Array[String] bam_addresses
        Int target_coverage
        String remote_dir
    }
    parameter_meta {
        remote_dir: "Root directory in the remote bucket. Every sample is stored in a subdirectory."
        sample_id_clone: "ID to assign to the clone"
    }
    
    call HPRCDownloadClonesImpl {
        input:
            sample_id_clone = sample_id_clone,
            bam_addresses = bam_addresses,
            target_coverage = target_coverage,
            remote_dir = remote_dir
    }
    
    output {
    }
}


task HPRCDownloadClonesImpl {
    input {
        String sample_id_clone
        Array[String] bam_addresses
        Int target_coverage
        String remote_dir
    }
    parameter_meta {
    }
    
    Int disk_size_gb = (3*target_coverage)*8 + 512
    Int mem_gb = 128
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
        
        # 1. Randomizing the order of the BAMs
        LIST_FILE=~{write_lines(bam_addresses)}
        shuf ${LIST_FILE} > randomized.txt
        cat randomized.txt
        TARGET_N_BYTES=$(( ~{target_coverage} * 3000000000 * 2 ))
        touch tmp1.fastq
        while read ADDRESS; do
            SUCCESS="0"
            if [[ ${ADDRESS} == gs://* ]]; then
                SUCCESS=$(gsutil -u ${BILLING_PROJECT} cp ${ADDRESS} . && echo 1 || echo 0)
            else
                SUCCESS=$(wget ${ADDRESS} && echo 1 || echo 0)
            fi
            if [[ ${SUCCESS} -eq 0 ]]; then
                continue
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
        done < randomized.txt
        
        # 2. Randomizing the order of the reads
        cp ~{docker_dir}/*.class .
        ${TIME_COMMAND} java FlattenFastq tmp1.fastq "SEPARATOR" tmp2.txt
        rm -f tmp1.fastq
        export TMPDIR="./terashuf_tmp"; mkdir ${TMPDIR}
        export MEMORY=~{mem_gb}
        ulimit -n 100000
        ${TIME_COMMAND} ~{docker_dir}/terashuf/terashuf < tmp2.txt > tmp3.txt
        rm -f tmp2.txt; rm -rf ${TMPDIR}
        head -c ${TARGET_N_BYTES} tmp3.txt > tmp4.txt
        N_ROWS=$(wc -l < tmp4.txt)
        rm -f tmp4.txt
        N_ROWS=$(( ${N_ROWS}*2 ))
        head -n ${N_ROWS} tmp3.txt | sed 's/SEPARATOR/\n/g' | gzip -1 > ~{sample_id_clone}.fastq.gz
        rm -f tmp3.txt
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id_clone}.fastq.gz ~{remote_dir}/~{sample_id_clone}/ && echo 0 || echo 1)
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
        cpu: 32
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
