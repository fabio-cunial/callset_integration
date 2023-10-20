version 1.0


# 
#
workflow HGSVC2DownloadAssemblies {
    input {
        Array[String] sample_ids
        Array[File] addresses
        String remote_dir
    }
    parameter_meta {
        remote_dir: "Root directory in the remote bucket. Every sample is stored in a subdirectory."
    }
    
    scatter(i in range(length(sample_ids))) {
        call HGSVC2DownloadAssembliesImpl {
            input:
                sample_id = sample_ids[i],
                addresses = addresses[i],
                remote_dir = remote_dir
        }
    }
    
    output {
    }
}


task HGSVC2DownloadAssembliesImpl {
    input {
        String sample_id
        File addresses
        String remote_dir
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
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        while read ADDRESS; do
            ${TIME_COMMAND} wget ${ADDRESS} &
            ${TIME_COMMAND} wget ${ADDRESS}.fai &
        done < ~{addresses}
        wait
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp '*.fasta*' ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
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
        docker: "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/hgsvc2:fc_hgsvc2_setup"
        cpu: 2
        memory: "8 GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}
