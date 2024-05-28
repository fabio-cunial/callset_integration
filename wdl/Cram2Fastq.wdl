version 1.0


#
workflow Cram2Fastq {
    input {
        File cram_file
        File reference_fa
        File reference_fai
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
    }
    call Cram2FastqImpl {
        input:
            cram_file = cram_file,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb
    }
    output {
        File reads1 = Cram2FastqImpl.reads1
        File reads2 = Cram2FastqImpl.reads2
    }
}


task Cram2FastqImpl {
    input {
        File cram_file
        File reference_fa
        File reference_fai
        Int n_cpus
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 4*ceil(size(cram_file,"GB")) + ceil(size(reference_fa,"GB")) + 512
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} samtools collate --threads ${N_THREADS} -u -o collated.bam ~{cram_file} ./prefix
        ${TIME_COMMAND} samtools fastq --threads ${N_THREADS} --reference ~{reference_fa} -n -0 /dev/null -s /dev/null -1 reads1.fastq.gz -2 reads2.fastq.gz collated.bam 
    >>>
    
    output {
        File reads1 = work_dir + "/reads1.fastq.gz"
        File reads2 = work_dir + "/reads2.fastq.gz"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
