#!/bin/bash
#
#
INPUT_DIR="/mnt/disks/disk-1/dipcall_bug/windows/sweep"
TOOL_NAME="dipcall"
MIN_FRACTION="0.9"
WINDOWS_FILE="windows.txt"


set -euxo pipefail

cd ${INPUT_DIR}
ls -d */ > ${WINDOWS_FILE}

while read WINDOW; do
    VALUES=$( tail -n 1 ${INPUT_DIR}/${WINDOW}${TOOL_NAME}/edges.csv )
    n_non_ref_edges=$( echo ${VALUES} | cut -d , -f 4 )
    n_non_ref_edges_covered=$( echo ${VALUES} | cut -d , -f 5 )
    THRESHOLD=$( echo "scale=2; ${MIN_FRACTION}*${n_non_ref_edges}" | bc )
    THRESHOLD=$( printf %3.f ${THRESHOLD} )
    if [ ${n_non_ref_edges_covered} -lt ${THRESHOLD} ]; then
        echo "${WINDOW},${n_non_ref_edges},${n_non_ref_edges_covered}"
    fi
done < ${WINDOWS_FILE}
rm -f ${WINDOWS_FILE}
