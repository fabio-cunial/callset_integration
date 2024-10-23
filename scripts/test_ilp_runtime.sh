#!/bin/bash
#
INPUT_DIR="./output/run"
N_WINDOWS="10"
INPUT_DIR_SMALL="./input_small"
OUTPUT_DIR="./output_small"
TIME_COMMAND="/usr/bin/time --verbose"
HAPESTRY_COMMAND="/mnt/disks/disk-1/compression/sv_merge/build/solve_from_directory"
N_THREADS="32"

# Building windows
rm -rf ${INPUT_DIR_SMALL} list.txt
ls ${INPUT_DIR} | sort --random-sort | head -n ${N_WINDOWS} > list.txt
while read WINDOW; do
    mkdir -p ${INPUT_DIR_SMALL}/${WINDOW}
    cp ${INPUT_DIR}/${WINDOW}/* ${INPUT_DIR_SMALL}/${WINDOW}/
done < list.txt

# Comparing compressed/non-compressed
COMPRESS=""
rm -rf ${OUTPUT_DIR} /tmp/* log_uncompressed.txt
${TIME_COMMAND} ${HAPESTRY_COMMAND} ${COMPRESS} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS} &> log_uncompressed.txt

COMPRESS="--compress_transmap true"
rm -rf ${OUTPUT_DIR} /tmp/* log_compressed.txt
${TIME_COMMAND} ${HAPESTRY_COMMAND} ${COMPRESS} --input ${INPUT_DIR_SMALL} --output_dir ${OUTPUT_DIR} --solver scip --n_threads ${N_THREADS} &> log_compressed.txt
