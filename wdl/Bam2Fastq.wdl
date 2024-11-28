version 1.0


#
workflow Bam2Fastq {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
    }
    
    call Bam2FastqImpl {
        input:
            sample_id = sample_id,
            input_bam = input_bam,
            input_bai = input_bai,
            n_cores = n_cores,
            mem_gb = mem_gb
    }
    
    output {
        File out_fastq = Bam2FastqImpl.out_fastq
    }
}


task Bam2FastqImpl {
    input {
        String sample_id
        File input_bam
        File input_bai
        Int n_cores
        Int mem_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int disk_size_gb = 10*ceil(size(input_bam, "GB")) + 256
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        df -h
        
        ${TIME_COMMAND} samtools collate --threads ${N_THREADS} -u -o collated.bam ~{input_bam} ./prefix
        samtools fastq -@ ${N_THREADS} -n collated.bam | pigz --processes ${N_THREADS} --fast --to-stdout > ~{sample_id}.fastq.gz
    >>>
    
    output {
        File out_fastq = work_dir + "/" + sample_id + ".fastq.gz"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: n_cores
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
