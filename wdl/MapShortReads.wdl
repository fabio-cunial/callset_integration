version 1.0


#
workflow MapShortReads {
    input {
        String remote_dir
        File reads1_fastq_gz
        File reads2_fastq_gz
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    call MapImpl {
        input:
            remote_dir = remote_dir,
            reads1_fastq_gz = reads1_fastq_gz,
            reads2_fastq_gz = reads2_fastq_gz,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    output {
        File alignments_bam = MapImpl.alignments_bam
        File alignments_bai = MapImpl.alignments_bai
    }
}


task MapImpl {
    input {
        String remote_dir
        File reads1_fastq_gz
        File reads2_fastq_gz
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(reads1_fastq_gz, "GB")) + ceil(size(reads2_fastq_gz, "GB")) + 10*ceil(size(reference_fa, "GB")) + 512
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
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
        
        # Downloading indexes
        REFERENCE_FILE=$(basename ~{reference_fa})
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{remote_dir}/${REFERENCE_FILE}.'*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading indexes. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        mv ~{reference_fa} ./${REFERENCE_FILE}
        mv ~{reference_fai} ./${REFERENCE_FILE}.fai
        
        # Mapping
        ${TIME_COMMAND} ~{docker_dir}/bwa/bwa mem -t ${N_THREADS} -Y -o alignments.sam ${REFERENCE_FILE} ~{reads1_fastq_gz} ~{reads2_fastq_gz}
        ${TIME_COMMAND} samtools sort -@ ${N_THREADS} --output-fmt BAM -o alignments_sorted.bam alignments.sam
        rm -f alignments.sam
        samtools index -@ ${N_THREADS} alignments_sorted.bam
    >>>
    
    output {
        File alignments_bam = work_dir + "/alignments_sorted.bam"
        File alignments_bai = work_dir + "/alignments_sorted.bam.bai"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: 64
        memory: "128GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
