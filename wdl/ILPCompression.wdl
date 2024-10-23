version 1.0


# 
#
workflow ILPCompression {
    input {
        File non_sequence_data_tar_gz
        Int n_windows
    }
    parameter_meta {
    }
    
    call ILPCompressionImpl {
        input:
            non_sequence_data_tar_gz = non_sequence_data_tar_gz,
            n_windows = n_windows
    }
    
    output {
        File monitor_log = ILPCompressionImpl.monitor_log
    }
}


#
task ILPCompressionImpl {
    input {
        File non_sequence_data_tar_gz
        Int n_windows
    }
    parameter_meta {
    }
    
    String docker_dir = "/hapestry"
    String work_dir = "/cromwell_root/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        INPUT_DIR="./output/run"
        INPUT_DIR_SMALL="./input_sample"
        OUTPUT_DIR="./output_sample"
        HAPESTRY_COMMAND="~{docker_dir}/sv_merge/build/solve_from_directory"

        # Sampling windows
        mv ~{non_sequence_data_tar_gz} ./archive.tar.gz
        tar -xzf ./archive.tar.gz
        rm -f ./archive.tar.gz
        rm -rf ${INPUT_DIR_SMALL} list.txt
        ls ${INPUT_DIR} | sort --random-sort | head -n ~{n_windows} > list.txt
        while read WINDOW; do
            mkdir -p ${INPUT_DIR_SMALL}/${WINDOW}
            mv ${INPUT_DIR}/${WINDOW}/* ${INPUT_DIR_SMALL}/${WINDOW}/
        done < list.txt
        rm -rf ${INPUT_DIR}
        
        # Starting resource monitoring
        bash ~{docker_dir}/vm_local_monitoring_script.sh &> monitoring.log &
        MONITOR_JOB=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')
        
        # Comparing compressed/non-compressed
        rm -rf ${OUTPUT_DIR} /tmp/* uncompressed.log
        ${TIME_COMMAND} ${HAPESTRY_COMMAND} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS}
        rm -rf ${OUTPUT_DIR} /tmp/* compressed.log
        ${TIME_COMMAND} ${HAPESTRY_COMMAND} --compress_transmap true --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS}
        
        # Stopping resource monitoring
        kill ${MONITOR_JOB}
        tail -n 100 ${MONITOR_FILE}
    >>>

    output {
        File monitor_log = work_dir + "/monitoring.log"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 32
        memory: "64GB"
        disks: "local-disk 500 HDD"
        preemptible: 0
    }
}
