version 1.0


#
workflow BwaIndex {
    input {
        String remote_dir
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    call IndexImpl {
        input:
            remote_dir = remote_dir,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    output {
    }
}


task IndexImpl {
    input {
        String remote_dir
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    Int ram_size_gb = 20*ceil(size(reference_fa, "GB"))
    Int disk_size_gb = 40*ceil(size(reference_fa, "GB"))
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
        
        REFERENCE_FILE=$(basename ~{reference_fa})
        ${TIME_COMMAND} ~{docker_dir}/bwa/bwa index ~{reference_fa}
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{reference_fa}'.*' ~{remote_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading indexes. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: 1
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
