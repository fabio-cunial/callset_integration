#!/bin/bash
INPUT_DIR=$1

set -euxo pipefail

TOOL1="bcftools_merge_chr1_norm"
TOOL2="svmerger_chr1_norm"

find ${INPUT_DIR} -type d -mindepth 1 -maxdepth 1 -print > directories.txt
while read DIRECTORY; do
    NODES1=$(grep ^S ${DIRECTORY}/${TOOL1}/graph.gfa | wc -l)
    NODES2=$(grep ^S ${DIRECTORY}/${TOOL2}/graph.gfa | wc -l)
    if [ ${NODES2} -gt ${NODES1} ]; then
        echo "ERROR at ${DIRECTORY}: ${NODES2} > ${NODES1}"
    fi
done < directories.txt
