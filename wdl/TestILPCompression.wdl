version 1.0


#
workflow TestILPCompression {
    input {
        File non_sequence_data_tar_gz
        Float d_weight = 4
        Float n_weight = 1
        Int max_timeout_minutes = 60
        Boolean run_uncompressed = true
        String gurobi_license = "null"
        String docker = "fcunial/hapestry:compression"
    }
    parameter_meta {
    }
    
    call TestILPCompressionImpl {
        input:
            non_sequence_data_tar_gz = non_sequence_data_tar_gz,
            max_timeout_minutes = max_timeout_minutes,
            d_weight = d_weight,
            n_weight = n_weight,
            run_uncompressed = run_uncompressed,
            gurobi_license = gurobi_license,
            docker = docker
    }
    
    output {
    }
}


task TestILPCompressionImpl {
    input {
        File non_sequence_data_tar_gz
        Float d_weight
        Float n_weight
        Int max_timeout_minutes
        Boolean run_uncompressed
        String gurobi_license
        String docker
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/hapestry"
    String docker_dir = "/hapestry"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        #N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        N_THREADS="1"
        HAPESTRY_COMMAND="~{docker_dir}/sv_merge/build/solve_from_directory"
        INPUT_DIR="./output/run"
        INPUT_DIR_SMALL="./input_small"
        OUTPUT_DIR="./output_small"
        IDENTICAL_LOG="./identical.log"
        
        # If a license is provided, we just place it in the conventional
        # location to be found automatically
        if [ ~{gurobi_license} != "null" ]; then
            gsutil -m cp ~{gurobi_license} /opt/gurobi/
            SOLVER_FLAG="--solver gurobi"
        else
            SOLVER_FLAG="--solver scip"
        fi
        
        tar -xzf ~{non_sequence_data_tar_gz}
        rm -rf ${INPUT_DIR_SMALL} list.txt ${IDENTICAL_LOG}
        ls ${INPUT_DIR} > list.txt && echo 0 || echo 1
        while read WINDOW; do
            # Building windows
            mkdir -p ${INPUT_DIR_SMALL}/${WINDOW}
            cp ${INPUT_DIR}/${WINDOW}/* ${INPUT_DIR_SMALL}/${WINDOW}/
            rm -f ${INPUT_DIR_SMALL}/${WINDOW}/solution.csv

            # Compressed
            rm -rf ${OUTPUT_DIR} log_compressed.txt && echo 0 || echo 1
            timeout ~{max_timeout_minutes}m ${HAPESTRY_COMMAND} --compress_transmap true --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} ${SOLVER_FLAG} --n_threads ${N_THREADS} --d_weight ~{d_weight} --n_weight ~{n_weight} &> log_compressed.txt || echo "0"
            if [ -e ${OUTPUT_DIR}/${WINDOW}/solution.csv ]; then
                mv ${OUTPUT_DIR}/${WINDOW}/solution.csv ${INPUT_DIR_SMALL}/${WINDOW}/solution_compressed.csv
            fi
            if [ -e log_compressed.txt ]; then
                cat log_compressed.txt
            fi
            echo "---- COMPRESSED log.csv:"
            if [ -e ${OUTPUT_DIR}/${WINDOW}/log.csv ]; then
                cat ${OUTPUT_DIR}/${WINDOW}/log.csv
            fi
            
            if ~{run_uncompressed} ; then
                # Uncompressed
                rm -rf ${OUTPUT_DIR} log_uncompressed.txt && echo 0 || echo 1
                timeout ~{max_timeout_minutes}m ${HAPESTRY_COMMAND} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} ${SOLVER_FLAG} --n_threads ${N_THREADS} --d_weight ~{d_weight} --n_weight ~{n_weight} &> log_uncompressed.txt || echo "0"
                if [ -e ${OUTPUT_DIR}/${WINDOW}/solution.csv ]; then
                    mv ${OUTPUT_DIR}/${WINDOW}/solution.csv ${INPUT_DIR_SMALL}/${WINDOW}/solution_uncompressed.csv
                fi
                echo "---- UNCOMPRESSED log.csv:"
                if [ -e ${OUTPUT_DIR}/${WINDOW}/log.csv ]; then
                    cat ${OUTPUT_DIR}/${WINDOW}/log.csv
                fi

                if [ ! -e ${INPUT_DIR_SMALL}/${WINDOW}/solution_uncompressed.csv -o ! -e ${INPUT_DIR_SMALL}/${WINDOW}/solution_compressed.csv ]; then
                    # Next iteration
                    rm -rf ${INPUT_DIR_SMALL}/${WINDOW}
                    continue
                fi

                # Checking if objective values are identical
                cut -d ',' -f 1,3 ${INPUT_DIR_SMALL}/${WINDOW}/solution_uncompressed.csv | sort | uniq > uncompressed_solution.txt
                cut -d ',' -f 1,3 ${INPUT_DIR_SMALL}/${WINDOW}/solution_compressed.csv | sort | uniq > compressed_solution.txt
                diff --brief uncompressed_solution.txt compressed_solution.txt &> test.txt || echo 0
                TEST=$(wc -l < test.txt)
                if [ ${TEST} -gt 0 ]; then
                    FIELD_1="0"
                else
                    FIELD_1="1"
                fi
                grep "Optimal" log_uncompressed.txt > uncompressed_optimum.txt
                grep "Optimal" log_compressed.txt > compressed_optimum.txt
                diff --brief uncompressed_optimum.txt compressed_optimum.txt &> test.txt || echo 0
                TEST=$(wc -l < test.txt)
                if [ ${TEST} -gt 0 ]; then
                    FIELD_2="0"
                    diff uncompressed_optimum.txt compressed_optimum.txt || echo 0
                else
                    FIELD_2="1"
                fi
                grep "Objective" log_uncompressed.txt > uncompressed_objective.txt
                grep "Objective" log_compressed.txt > compressed_objective.txt
                diff --brief uncompressed_objective.txt compressed_objective.txt &> test.txt || echo 0
                TEST=$(wc -l < test.txt)
                if [ ${TEST} -gt 0 ]; then
                    FIELD_3="0"
                    diff uncompressed_objective.txt compressed_objective.txt || echo 0
                else
                    FIELD_3="1"
                fi
                grep "d_min" log_uncompressed.txt
                grep "d_min" log_compressed.txt
                grep "n_max" log_uncompressed.txt
                grep "n_max" log_compressed.txt
                if [ ${FIELD_2} -eq 0 -o ${FIELD_3} -eq 0 ]; then
                    echo "ERROR: compressed and uncompressed objectives differ in window ${WINDOW}"
                    exit 1
                fi
                echo "${FIELD_1},${FIELD_2},${FIELD_3}" >> ${IDENTICAL_LOG}
            fi

            # Next iteration
            rm -rf ${INPUT_DIR_SMALL}/${WINDOW}
        done < list.txt
    >>>

    output {
    }
    runtime {
        docker: docker
        cpu: 1
        memory: "64GB"
        disks: "local-disk 500 HDD"
        preemptible: 0
    }
}