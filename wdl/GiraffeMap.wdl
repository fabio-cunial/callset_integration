version 1.0


#
workflow GiraffeMap {
    input {
        String remote_dir
        String sample_id
        File reads1_fastq_gz
        File reads2_fastq_gz
        Int n_cpus = 16
        Int ram_size_gb = 32
        Int disk_size_gb
    }
    parameter_meta {
        remote_dir: "The GBZ, MIN and DIST graph indexes are assumed to be in directory `remote_dir/sample_id`."
    }
    call GiraffeMapImpl {
        input:
            remote_dir = remote_dir,
            sample_id = sample_id,
            reads1_fastq_gz = reads1_fastq_gz,
            reads2_fastq_gz = reads2_fastq_gz,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File alignments_gam = GiraffeMapImpl.alignments_gam
        File stats = GiraffeMapImpl.stats
    }
}


# COMMAND    | TIME | CORES | RAM
# vg giraffe | 6h   |   16  | 50G
# vg stats   |  ?   |    ?  |   ?     CRASHED
#
task GiraffeMapImpl {
    input {
        String remote_dir
        String sample_id
        File reads1_fastq_gz
        File reads2_fastq_gz
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
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
        VG_COMMAND="~{docker_dir}/vg"
        
        
        while : ; do
            TEST=$(gsutil -m cp ~{remote_dir}/~{sample_id}/~{sample_id}.gbz ~{remote_dir}/~{sample_id}/~{sample_id}.min ~{remote_dir}/~{sample_id}/~{sample_id}.dist . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading indexes. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        ${TIME_COMMAND} ${VG_COMMAND} giraffe --threads ${N_THREADS} --progress --output-format gam --gbz-name ~{sample_id}.gbz --minimizer-name ~{sample_id}.min --dist-name ~{sample_id}.dist --fastq-in ~{reads1_fastq_gz} --fastq-in ~{reads2_fastq_gz} > ~{sample_id}.gam
        # The following command crashes on vg 1.58.0:
        #${TIME_COMMAND} ${VG_COMMAND} stats --threads ${N_THREADS} --alignments ~{sample_id}.gam > stats.txt
        touch stats.txt
        df -h
    >>>
    
    output {
        File alignments_gam = work_dir + "/" + sample_id + ".gam"
        File stats = work_dir + "/stats.txt"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
