version 1.0


# 
#
workflow ILPCompression {
    input {
        String remote_run_dir
        Int min_minutes = 1
        Int max_minutes = 14
        Int n_windows = 50
    }
    parameter_meta {
    }
    
    call ILPCompressionImpl {
        input:
            remote_run_dir = remote_run_dir,
            min_minutes = min_minutes,
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
        String remote_run_dir
        Int min_minutes
        Int max_minutes
        Int n_windows
    }
    parameter_meta {
        remote_run_dir: "The parent directory of all the `shard-X` directories."
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
        MAX_CHUNK_ID="128"
        
        
        # Returns 1 iff the window is sampled, i.e. if all the ILPs took
        # $<=max_minutes$ and at least one ILP took $>=min_minutes$.
        #
        function sampleWindow() {
            local WINDOW=$1
    
            LOG_FILE="${WINDOW}/log.csv"
            if [ ! -e ${LOG_FILE} ]; then
                echo "0"
                return 0
            fi
            
            MAX_MINS="0"
            
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
            if [ ${MINUTES_F} -gt ${MAX_MINS} ]; then
                MAX_MINS=${MINUTES_F}
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
            if [ ${MINUTES_D} -gt ${MAX_MINS} ]; then
                MAX_MINS=${MINUTES_D}
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
            if [ ${MINUTES_ND} -gt ${MAX_MINS} ]; then
                MAX_MINS=${MINUTES_ND}
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
            if [ ${MINUTES_DN} -gt ${MAX_MINS} ]; then
                MAX_MINS=${MINUTES_DN}
            fi
            
            if [ ${MAX_MINS} -ge ~{min_minutes} ]; then
                echo "1"
            else
                echo "0"
            fi
        }
        
        
        # Sampling $n_windows$ that succeeded and that took $<=max_minutes$ in
        # every ILP.
        CHUNK_ID="0"; N_SAMPLED_WINDOWS="0"
        while [ ${N_SAMPLED_WINDOWS} -lt ~{n_windows} ]; do
            TEST=$(gsutil ls ~{remote_run_dir}/shard-${CHUNK_ID}/output/non_sequence_data.tar.gz)
            if [ -n ${TEST} ]; then
                gsutil -m cp ~{remote_run_dir}/shard-${CHUNK_ID}/output/non_sequence_data.tar.gz .
            else
                TEST=$(gsutil ls ~{remote_run_dir}/shard-${CHUNK_ID}/attempt-2/output/non_sequence_data.tar.gz)
                if [ -n ${TEST} ]; then
                    gsutil -m cp ~{remote_run_dir}/shard-${CHUNK_ID}/attempt-2/output/non_sequence_data.tar.gz .
                else 
                    CHUNK_ID=$(( ${CHUNK_ID} + 1 ))
                    if [ ${CHUNK_ID} -gt ${MAX_CHUNK_ID} ]; then
                        break
                    else
                        continue
                    fi
                fi
            fi
            tar -xzf ./non_sequence_data.tar.gz
            rm -f non_sequence_data.tar.gz
            ls ${INPUT_DIR} | sort --random-sort > list.txt && echo 0 || echo 1
            while read WINDOW; do
                cat ${INPUT_DIR}/${WINDOW}/log.csv
                SUCCESS=$(sampleWindow ${INPUT_DIR}/${WINDOW})
                if [ ${SUCCESS} -eq 1 ]; then
                    mkdir -p ${INPUT_DIR_SMALL}/${WINDOW}
                    mv ${INPUT_DIR}/${WINDOW}/* ${INPUT_DIR_SMALL}/${WINDOW}/
                    N_SAMPLED_WINDOWS=$(( ${N_SAMPLED_WINDOWS} + 1 ))
                    if [ ${N_SAMPLED_WINDOWS} -eq ~{n_windows} ]; then
                        break
                    fi
                fi
            done < list.txt
            rm -rf ${INPUT_DIR}
            CHUNK_ID=$(( ${CHUNK_ID} + 1 ))
            if [ ${CHUNK_ID} -gt ${MAX_CHUNK_ID} ]; then
                break
            fi
        done
        
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
