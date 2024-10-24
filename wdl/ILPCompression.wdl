version 1.0


# 
#
workflow ILPCompression {
    input {
        File non_sequence_data_tar_gz
        Int max_minutes = 14
        Int n_windows
    }
    parameter_meta {
    }
    
    call ILPCompressionImpl {
        input:
            non_sequence_data_tar_gz = non_sequence_data_tar_gz,
            max_minutes = max_minutes,
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
        Int max_minutes
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
        
        
        # Returns 1 iff the window is sampled.
        #
        function sampleWindow() {
            local WINDOW=$1
    
            LOG_FILE="${WINDOW}/log.csv"
            if [ ! -e ${LOG_FILE} ]; then
                echo "0"
                return 0
            fi
            cat ${LOG_FILE}
    
            grep 'feasibility,' ${LOG_FILE} > tmp.txt || echo ""
            N_LINES=$(wc -l < tmp.txt)
            if [ ${N_LINES} -eq 0 ]; then
                echo "0"
                return 0
            fi
            HOURS_F=$( cat tmp.txt | cut -d , -f 2 )
            MINUTES_F=$( cat tmp.txt | cut -d , -f 3 )
            SUCCESS_F=$( cat tmp.txt | cut -d , -f 6 )
            if [ ${SUCCESS_F} -eq 0 -o ${HOURS_F} -ge 1 -o ${MINUTES_F} -gt ~{max_minutes} ]; then
                echo "0"
                return 0
            fi
    
            grep 'optimize_d,' ${LOG_FILE} > tmp.txt || echo ""
            N_LINES=$(wc -l < tmp.txt)
            if [ ${N_LINES} -eq 0 ]; then
                echo "0"
                return 0
            fi
            HOURS_D=$( cat tmp.txt | cut -d , -f 2 )
            MINUTES_D=$( cat tmp.txt | cut -d , -f 3 )
            SUCCESS_D=$( cat tmp.txt | cut -d , -f 6 )
            if [ ${SUCCESS_F} -eq 0 -o ${HOURS_F} -ge 1 -o ${MINUTES_F} -gt ~{max_minutes} ]; then
                echo "0"
                return 0
            fi
    
            grep 'optimize_n_given_d,' ${LOG_FILE} > tmp.txt || echo ""
            N_LINES=$(wc -l < tmp.txt)
            if [ ${N_LINES} -eq 0 ]; then
                echo "0"
                return 0
            fi
            HOURS_ND=$( cat tmp.txt | cut -d , -f 2 )
            MINUTES_ND=$( cat tmp.txt | cut -d , -f 3 )
            SUCCESS_ND=$( cat tmp.txt | cut -d , -f 6 )
            if [ ${SUCCESS_F} -eq 0 -o ${HOURS_F} -ge 1 -o ${MINUTES_F} -gt ~{max_minutes} ]; then
                echo "0"
                return 0
            fi
    
            grep 'optimize_d_plus_n,' ${LOG_FILE} > tmp.txt || echo ""
            N_LINES=$(wc -l < tmp.txt)
            if [ ${N_LINES} -eq 0 ]; then
                echo "0"
                return 0
            fi
            HOURS_DN=$( cat tmp.txt | cut -d , -f 2 )
            MINUTES_DN=$( cat tmp.txt | cut -d , -f 3 )
            SUCCESS_DN=$( cat tmp.txt | cut -d , -f 6 )
            if [ ${SUCCESS_F} -eq 0 -o ${HOURS_F} -ge 1 -o ${MINUTES_F} -gt ~{max_minutes} ]; then
                echo "0"
                return 0
            fi
            
            echo "1"
        }
        
        
        # Sampling $n_windows$ that succeeded and that took $<=max_minutes$ in
        # every ILP.
        mv ~{non_sequence_data_tar_gz} ./archive.tar.gz
        tar -xzf ./archive.tar.gz
        rm -f ./archive.tar.gz
        ls ${INPUT_DIR} | sort --random-sort > list.txt && echo 0 || echo 1
        N_WINDOWS="0"
        while read WINDOW; do
            SUCCESS=$(sampleWindow ${INPUT_DIR}/${WINDOW})
            if [ ${SUCCESS} -eq 1 ]; then
                mkdir -p ${INPUT_DIR_SMALL}/${WINDOW}
                mv ${INPUT_DIR}/${WINDOW}/* ${INPUT_DIR_SMALL}/${WINDOW}/
                N_WINDOWS=$(( ${N_WINDOWS} + 1 ))
                if [ ${N_WINDOWS} -eq ~{n_windows} ]; then
                    break
                fi
            fi
        done < list.txt
        rm -rf ${INPUT_DIR}
        
        # Starting resource monitoring
        export MONITOR_MOUNT_POINT=~{work_dir}
        bash ~{docker_dir}/vm_local_monitoring_script.sh &> monitoring.log &
        MONITOR_JOB=$(ps -aux | grep -F 'vm_local_monitoring_script.sh' | head -1 | awk '{print $2}')
        
        # Comparing compressed/non-compressed
        rm -rf ${OUTPUT_DIR} /tmp/*
        ${TIME_COMMAND} ${HAPESTRY_COMMAND} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS}
        rm -rf ${OUTPUT_DIR} /tmp/*
        ${TIME_COMMAND} ${HAPESTRY_COMMAND} --compress_transmap true --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS}
        
        # Stopping resource monitoring
        kill ${MONITOR_JOB}
        tail -n 100 monitoring.log
    >>>

    output {
        File monitor_log = work_dir + "/monitoring.log"
    }
    runtime {
        docker: "fcunial/hapestry:compression"
        cpu: 32
        memory: "64GB"
        disks: "local-disk 500 HDD"
        preemptible: 0
    }
}
