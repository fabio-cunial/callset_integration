#!/bin/bash
#
INPUT_DIR="./problematic_window"   #"./problematic_window"  #"./output/run"
N_WINDOWS="200"
INPUT_DIR_SMALL="./input_small"
OUTPUT_DIR="./output_small"
TIME_COMMAND="/usr/bin/time --verbose"
HAPESTRY_COMMAND="/mnt/disks/disk-1/compression/sv_merge/build/solve_from_directory"
N_THREADS="1"
IDENTICAL_LOG="./identical.log"
MAX_MINUTES="2"
MAX_TIMEOUT="2"

set -euxo pipefail

rm -rf ${INPUT_DIR_SMALL} list.txt ${IDENTICAL_LOG}
ls ${INPUT_DIR} | sort --random-sort | head -n ${N_WINDOWS} > list.txt && echo 0 || echo 1
while read WINDOW; do
    # Building windows
    mkdir -p ${INPUT_DIR_SMALL}/${WINDOW}
    cp ${INPUT_DIR}/${WINDOW}/* ${INPUT_DIR_SMALL}/${WINDOW}/ 
    rm -f ${INPUT_DIR_SMALL}/${WINDOW}/solution.csv

    # Uncompressed SCIP
    rm -rf ${OUTPUT_DIR} log_scip.txt && echo 0 || echo 1
    timeout ${MAX_TIMEOUT}m ${HAPESTRY_COMMAND} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS} &> log_scip.txt || echo "0"
    if [ -e ${OUTPUT_DIR}/${WINDOW}/solution.csv ]; then 
        mv ${OUTPUT_DIR}/${WINDOW}/solution.csv ${INPUT_DIR_SMALL}/${WINDOW}/solution_scip.csv
    fi
    echo "---- UNCOMPRESSED SCIP log.csv:"
    cat ${OUTPUT_DIR}/${WINDOW}/log.csv
    
    # Uncompressed GUROBI
    rm -rf ${OUTPUT_DIR} log_gurobi.txt && echo 0 || echo 1
    timeout ${MAX_TIMEOUT}m ${HAPESTRY_COMMAND} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver gurobi --n_threads ${N_THREADS} &> log_gurobi.txt || echo "0"
    if [ -e ${OUTPUT_DIR}/${WINDOW}/solution.csv ]; then 
        mv ${OUTPUT_DIR}/${WINDOW}/solution.csv ${INPUT_DIR_SMALL}/${WINDOW}/solution_gurobi.csv
    fi
    echo "---- UNCOMPRESSED GUROBI log.csv:"
    cat ${OUTPUT_DIR}/${WINDOW}/log.csv

    # Checking if objective values are identical
    cut -d ',' -f 1,3 ${INPUT_DIR_SMALL}/${WINDOW}/solution_scip.csv | sort | uniq > scip_solution.txt
    cut -d ',' -f 1,3 ${INPUT_DIR_SMALL}/${WINDOW}/solution_gurobi.csv | sort | uniq > gurobi_solution.txt
    diff --brief scip_solution.txt gurobi_solution.txt &> test.txt || echo 0
    TEST=$(wc -l < test.txt)
    if [ ${TEST} -gt 0 ]; then
        FIELD_1="0"
    else
        FIELD_1="1"
    fi
    grep "Optimal" log_scip.txt > scip_optimum.txt
    grep "Optimal" log_gurobi.txt > gurobi_optimum.txt
    diff --brief scip_optimum.txt gurobi_optimum.txt &> test.txt || echo 0
    TEST=$(wc -l < test.txt)
    if [ ${TEST} -gt 0 ]; then
        FIELD_2="0"
    else
        FIELD_2="1"
    fi
    grep "Objective" log_scip.txt > scip_objective.txt
    grep "Objective" log_gurobi.txt > gurobi_objective.txt
    diff --brief scip_objective.txt gurobi_objective.txt &> test.txt || echo 0
    TEST=$(wc -l < test.txt)
    if [ ${TEST} -gt 0 ]; then
        FIELD_3="0"
    else
        FIELD_3="1"
    fi
    grep "d_min" log_scip.txt
    grep "d_min" log_gurobi.txt
    grep "n_max" log_scip.txt
    grep "n_max" log_gurobi.txt
    if [ ${FIELD_2} -eq 0 -o ${FIELD_3} -eq 0 ]; then
        echo "ERROR: SCIP and GUROBI objectives differ in window ${WINDOW}"
        cat scip_objective.txt gurobi_objective.txt
        exit 1
    fi
    echo "${FIELD_1},${FIELD_2},${FIELD_3}" >> ${IDENTICAL_LOG}
    
    # Next iteration
    rm -rf ${INPUT_DIR_SMALL}/${WINDOW}
done < list.txt
