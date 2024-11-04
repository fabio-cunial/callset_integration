#!/bin/bash
#
INPUT_DIR="./output/run"
N_WINDOWS="200"
INPUT_DIR_SMALL="./input_small"
OUTPUT_DIR="./output_small"
TIME_COMMAND="/usr/bin/time --verbose"
HAPESTRY_COMMAND="/mnt/disks/disk-1/compression/sv_merge/build/solve_from_directory"
N_THREADS="1"
IDENTICAL_LOG="./identical.log"

set -euxo pipefail

rm -rf ${INPUT_DIR_SMALL} list.txt ${IDENTICAL_LOG}
ls ${INPUT_DIR} | sort --random-sort | head -n ${N_WINDOWS} > list.txt && echo 0 || echo 1
while read WINDOW; do
    # Building windows
    mkdir -p ${INPUT_DIR_SMALL}/${WINDOW}
    cp ${INPUT_DIR}/${WINDOW}/* ${INPUT_DIR_SMALL}/${WINDOW}/ 
    rm -f ${INPUT_DIR_SMALL}/${WINDOW}/solution.csv
    
    # Uncommpressed
    rm -rf ${OUTPUT_DIR} log_uncompressed.txt && echo 0 || echo 1
    ${TIME_COMMAND} ${HAPESTRY_COMMAND} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS} &> log_uncompressed.txt
    if [ -e ${OUTPUT_DIR}/${WINDOW}/solution.csv ]; then 
        mv ${OUTPUT_DIR}/${WINDOW}/solution.csv ${INPUT_DIR_SMALL}/${WINDOW}/solution_uncompressed.csv
    fi

    # Compressed
    rm -rf ${OUTPUT_DIR} log_compressed.txt && echo 0 || echo 1
    ${TIME_COMMAND} ${HAPESTRY_COMMAND} --compress_transmap true --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS} &> log_compressed.txt
    if [ -e ${OUTPUT_DIR}/${WINDOW}/solution.csv ]; then 
        mv ${OUTPUT_DIR}/${WINDOW}/solution.csv ${INPUT_DIR_SMALL}/${WINDOW}/solution_compressed.csv
    fi
    
    # Checking if solutions are identical
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
    else
        FIELD_2="1"
    fi
    grep "Objective" log_uncompressed.txt > uncompressed_objective.txt
    grep "Objective" log_compressed.txt > compressed_objective.txt
    diff --brief uncompressed_objective.txt compressed_objective.txt &> test.txt || echo 0
    TEST=$(wc -l < test.txt)
    if [ ${TEST} -gt 0 ]; then
        FIELD_3="0"
    else
        FIELD_3="1"
    fi
#    grep "d_min" log_uncompressed.txt
#    grep "d_min" log_compressed.txt
#    grep "n_max" log_uncompressed.txt
#    grep "n_max" log_compressed.txt
    if [ ${FIELD_1} -eq 0 -a ${FIELD_3} -eq 0 ]; then
        echo "ERROR: compressed and uncompressed differ in window ${WINDOW}"
        exit 1
    fi
    echo "${FIELD_1},${FIELD_2},${FIELD_3}" >> ${IDENTICAL_LOG}
    
    # Next iteration
    rm -rf ${INPUT_DIR_SMALL}/${WINDOW}
done < list.txt
