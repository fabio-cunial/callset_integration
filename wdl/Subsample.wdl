version 1.0


# Downloads BAM files up to a given coverage, randomizing both the order in
# which the BAMs are downloaded, and the order of the reads in the
# concatenation of all downloaded BAMs.
#
# Remark: the output FASTQ file for coverage $i+1$ is a superset of the output
# FASTQ file for coverage $i$.
#
# Performance on HPRC's HG002. Requested coverages: 8,16,32.
# 16 physical cores, 128 GB RAM, 2 TB HDD. 
# Total time: 3h 30m
# Total cost: ????
#
# STEP                  TIME            CPU            RAM
# samtools fastq         1 m            250 %         16 MB
# seqkit stats          13 s            100 %         20 MB
# seqkit scat           40 m            150 %          4 GB
# split                 25 m             10 %          2 MB
# FlattenFastq           9 m             10 %          4 GB
# terashuf              27 m             23 %        116 GB
# UnflattenFastq        30 m             30 %          4 GB
# pigz                  11 m            250 %         15 MB
#
workflow Subsample {
    input {
        String sample_id
        Array[String] bam_addresses
        String coverages
        String remote_dir
        String billing_project = "broad-firecloud-dsde-methods"
        Int haploid_genome_length_gb = 3
        Int n_cores = 16
        Int mem_gb = 128
        Int disk_size_gb = 500
    }
    parameter_meta {
        bam_addresses: "Can be .bam, .fastq, .fastq.gz"
        coverages: "Comma-separated. Example: 8,16,32"
        remote_dir: "Output directory in a remote bucket"
        n_cores: ">=max{4, 2*n_coverages}"
    }
    
    call SubsampleImpl {
        input:
            sample_id = sample_id,
            bam_addresses = bam_addresses,
            coverages = coverages,
            remote_dir = remote_dir,
            billing_project = billing_project,
            haploid_genome_length_gb = haploid_genome_length_gb,
            n_cores = n_cores,
            mem_gb = mem_gb,
            disk_size_gb = disk_size_gb
    }
    
    output {
    }
}


task SubsampleImpl {
    input {
        String sample_id
        Array[String] bam_addresses
        String coverages
        String remote_dir
        String billing_project
        Int haploid_genome_length_gb
        Int n_cores
        Int mem_gb
        Int disk_size_gb
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
        RAM_PER_THREAD_BYTES=$(( (1000000000*( ~{mem_gb} -5)) / ${N_THREADS} ))
        HAPLOID_GENOME_LENGTH_GB=$(( ~{haploid_genome_length_gb} * 1000000000 ))
        df -h
        
        
        function flattenFastq() {
            local ID=$1
            
            ${TIME_COMMAND} java FlattenFastq ${ID} "SEPARATOR" ${ID}.flattened ${RAM_PER_THREAD_BYTES}
            rm -f ${ID}
        }
        
        
        MAX_COVERAGE=~{coverages}
        MAX_COVERAGE=${MAX_COVERAGE##*,}
        
        # 1. Randomizing the order of the BAMs and downloading enough of them
        LIST_FILE=~{write_lines(bam_addresses)}
        shuf ${LIST_FILE} > randomized.txt
        cat randomized.txt
        TARGET_N_CHARS=$(( ${MAX_COVERAGE} * ${HAPLOID_GENOME_LENGTH_GB} ))
        TOTAL_N_CHARS="0"; TOTAL_N_READS="0"; FILE_ID="0";
        while read ADDRESS; do
            FILE_ID=$(( ${FILE_ID} + 1 ))
            SUCCESS="0"
            if [[ ${ADDRESS} == gs://* ]]; then
                SUCCESS=$(gsutil -u ~{billing_project} cp ${ADDRESS} . && echo 1 || echo 0)
            else
                SUCCESS=$(wget ${ADDRESS} && echo 1 || echo 0)
            fi
            if [[ ${SUCCESS} -eq 0 ]]; then
                continue
            fi
            FILE_NAME=$(basename ${ADDRESS})
            if [[ ${FILE_NAME} == *.bam ]]; then
                ${TIME_COMMAND} samtools fastq -@ ${N_THREADS} -n ${FILE_NAME} > ${FILE_ID}.fastq 
            elif [[ ${FILE_NAME} == *.fastq.gz ]]; then
                ${TIME_COMMAND} gunzip -c ${FILE_NAME} > ${FILE_ID}.fastq
            elif [[ ${FILE_NAME} == *.fastq ]]; then
                mv ${FILE_NAME} ${FILE_ID}.fastq
            fi
            rm -f ${FILE_NAME}
            ${TIME_COMMAND} ~{docker_dir}/seqkit stats --threads ${N_THREADS} --tabular ${FILE_ID}.fastq > stats.txt
            N_CHARS=$(cut -f 5 stats.txt | tail -n 1)
            TOTAL_N_CHARS=$(( ${TOTAL_N_CHARS} + ${N_CHARS} ))
            N_READS=$(cut -f 4 stats.txt | tail -n 1)
            TOTAL_N_READS=$(( ${TOTAL_N_READS} + ${N_READS} ))
            if [[ ${TOTAL_N_CHARS} -gt ${TARGET_N_CHARS} ]]; then
                break
            fi
            df -h
        done < randomized.txt
        mkdir ./fastqs
        mv *.fastq ./fastqs
        ${TIME_COMMAND} ~{docker_dir}/seqkit scat --threads ${N_THREADS} -f --out-format fastq ./fastqs > tmp1.fastq
        rm -rf ./fastqs
        df -h
        
        # 2. Randomizing the order of the reads
        LINES_PER_CHUNK=$((  ( ${TOTAL_N_READS} / ${N_THREADS} )*4  ))
        ${TIME_COMMAND} split -d -a 2 -l ${LINES_PER_CHUNK} tmp1.fastq
        rm -f tmp1.fastq
        cp ~{docker_dir}/*.class .
        for ID in $(seq -f '%02g' 0 $(( ${N_THREADS} + 1 )) ); do
            if [[ -e x${ID} ]]; then
                flattenFastq x${ID} &
            fi
        done
        wait
        cat *.flattened > tmp2.txt
        rm -f *.flattened
        export TMPDIR="./terashuf_tmp"; mkdir ${TMPDIR}
        export MEMORY=$(( ~{mem_gb} - 5 ))
        ulimit -n 100000
        ${TIME_COMMAND} ~{docker_dir}/terashuf/terashuf < tmp2.txt > tmp3.txt
        rm -f tmp2.txt; rm -rf ${TMPDIR}
        ${TIME_COMMAND} java UnflattenFastq tmp3.txt "SEPARATOR" ${HAPLOID_GENOME_LENGTH_GB} ~{coverages} tmp4.fastq 1000000000 > lastLines.txt
        rm -f tmp3.txt
        while read ROW; do
            if [[ -z ${ROW} ]]; then
                break
            fi
            COVERAGE=$(echo ${ROW} | cut -d , -f 1)
            LAST_LINE=$(echo ${ROW} | cut -d , -f 2)
            head -n ${LAST_LINE} tmp4.fastq | ${TIME_COMMAND} pigz --processes ${N_THREADS} --fast --to-stdout > ~{sample_id}_${COVERAGE}.fastq.gz
        done < lastLines.txt
        wait
        rm -f tmp4.fastq
        
        # 3. Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}_'*.fastq.gz' ~{remote_dir} && echo 0 || echo 1)
            if [[ ${TEST} -eq 1 ]]; then
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
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
