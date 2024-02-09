#!/bin/bash
INPUT_DIR=$1

set -euo pipefail

TOOL1="bcftools_merge_chr1_norm"
TOOL2="truvari_chr1_norm"

find ${INPUT_DIR} -type d -mindepth 1 -maxdepth 1 -print > directories.txt
while read DIRECTORY; do
    ALIGNMENTS1=$(wc -l < ${DIRECTORY}/${TOOL1}/alignments.gaf)
    ALIGNMENTS2=$(wc -l < ${DIRECTORY}/${TOOL2}/alignments.gaf)
    if [ ${ALIGNMENTS1} -gt ${ALIGNMENTS2} ]; then
        echo "ERROR at ${DIRECTORY}: ${ALIGNMENTS1} > ${ALIGNMENTS2}"
    fi
done < directories.txt
