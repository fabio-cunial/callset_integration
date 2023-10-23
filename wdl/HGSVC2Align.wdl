version 1.0


# Similar to an AoU production run. See:
# https://github.com/broadinstitute/long-read-pipelines/blob/8c4e9cdda7944d2e2ad36dce2bedd2be5e030fa0/wdl/tasks/PBUtils.wdl#L771
#
workflow HGSVC2Align {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        String remote_dir
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
        remote_dir: "Root directory in the remote bucket. Every sample is stored in a subdirectory."
    }
    
    call HGSVC2AlignImpl {
        input:
            sample_id = sample_id,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            reads_fastq_gz = reads_fastq_gz,
            remote_dir = remote_dir,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb
    }
    
    output {
    }
}


task HGSVC2AlignImpl {
    input {
        String sample_id
        File reference_fa
        File reference_fai
        File reads_fastq_gz
        String remote_dir
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(reads_fastq_gz, "GB")) + ceil(size(reference_fa, "GB")) + 100
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
        
        # Remark: pbmm2 automatically uses all cores.
        ${TIME_COMMAND} pbmm2 align --preset CCS --sort --sample ~{sample_id} ~{reference_fa} ~{reads_fastq_gz} out.bam
        ${TIME_COMMAND} samtools calmd -@ ${N_THREADS} --no-PG -b out.bam ~{reference_fa} > ~{sample_id}.bam
        ${TIME_COMMAND} samtools index -@ ${N_THREADS} ~{sample_id}.bam
        
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.'bam*' ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
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
