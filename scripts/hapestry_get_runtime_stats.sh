#!/bin/bash
#
WITH_SNPS_DIR=$1   # "/mnt/disks/disk-1/runtime_analysis/with_snps"
WITHOUT_SNPS_DIR=$2


LOG_CSV="log_csv.txt"
LOG_OPTIMIZER="log_optimizer.txt"
LOG_NPATHS="log_npaths.txt"
LOG_NEDGES="log_nedges.txt"

set -euo pipefail

rm -f ${WITH_SNPS_DIR}/${LOG_CSV} ${WITH_SNPS_DIR}/${LOG_OPTIMIZER} ${WITH_SNPS_DIR}/${LOG_NPATHS} ${WITH_SNPS_DIR}/${LOG_NEDGES}
rm -f ${WITHOUT_SNPS_DIR}/${LOG_CSV} ${WITHOUT_SNPS_DIR}/${LOG_OPTIMIZER} ${WITHOUT_SNPS_DIR}/${LOG_NPATHS} ${WITHOUT_SNPS_DIR}/${LOG_NEDGES}

# Basic consistency checks.
# Remark: a window in the without-SNPs file could have its boundaries changed in
# the with-SNPs file, since the with-SNPs file was built bu adding calls to the
# without-SNPs file.
ls ${WITH_SNPS_DIR}/output/run | sort > ${WITH_SNPS_DIR}/list.txt
N_WINDOWS_WITH_SNPS=$(wc -l < ${WITH_SNPS_DIR}/list.txt)
echo "WITH SNPS: ${N_WINDOWS_WITH_SNPS} windows"
ls ${WITHOUT_SNPS_DIR}/output/run | sort > ${WITHOUT_SNPS_DIR}/list.txt
N_WINDOWS_WITHOUT_SNPS=$(wc -l < ${WITHOUT_SNPS_DIR}/list.txt)
echo "WITHOUT SNPS: ${N_WINDOWS_WITHOUT_SNPS} windows"
if [ ${N_WINDOWS_WITHOUT_SNPS} -gt ${N_WINDOWS_WITH_SNPS} ]; then
    echo "Error: more windows without SNPs?!"
    exit
fi
comm -1 -2 ${WITH_SNPS_DIR}/list.txt ${WITHOUT_SNPS_DIR}/list.txt > shared_windows.txt
N_SHARED_WINDOWS=$(wc -l < shared_windows.txt)
echo "SHARED: ${N_SHARED_WINDOWS} windows"


function collectWindow() {
    local WINDOW=$1
    local LOG_CSV=$2
    local LOG_OPTIMIZER=$3
    local LOG_NPATHS=$4
    local LOG_NEDGES=$5
    
    # Reading log.csv
    LOG_FILE="${WINDOW}/log.csv"
    if [ -e ${LOG_FILE} ]; then
        HOURS=$( tail -n +2 ${LOG_FILE} | cut -d , -f 2 )
        MINUTES=$( tail -n +2 ${LOG_FILE} | cut -d , -f 3 )
        SECONDS=$( tail -n +2 ${LOG_FILE} | cut -d , -f 4 )
        SUCCESS=$( tail -n +2 ${LOG_FILE} | cut -d , -f 6 )
        TOTAL_SECONDS=$(( ${HOURS}*60*60 + ${MINUTES}*60 + ${SECONDS} ))
        if [ ${SUCCESS} -eq 1 ]; then
            echo "${TOTAL_SECONDS},1" >> ${LOG_CSV}
        else
            echo "${TOTAL_SECONDS},0" >> ${LOG_CSV}
        fi
    else
        echo "WARNING: log.csv missing from ${WINDOW}"
    fi
    
    # Reading log_optimizer.txt
    LOG_FILE="${WINDOW}/log_optimizer.txt"
    if [ -e ${LOG_FILE} ]; then
        TIME=$(tail -n 1 ${LOG_FILE} | cut -d ' ' -f 23)
        TIME=${TIME%,}
        echo ${TIME} >> ${LOG_OPTIMIZER}
    else
        echo "WARNING: log_optimizer.txt missing from ${WINDOW}"
    fi
    
    # Reading reads_to_paths.csv
    LOG_FILE="${WINDOW}/reads_to_paths.csv"
    if [ -e ${LOG_FILE} ]; then
        N_PATHS=$(cut -d , -f 4 ${LOG_FILE} | sort | uniq | wc -l)
        echo $(( ${N_PATHS} - 1 )) >> ${LOG_NPATHS}
        N_EDGES=$(wc -l < ${LOG_FILE})
        echo $(( ${N_EDGES} - 1 )) >> ${LOG_NEDGES}
    else
        echo "WARNING: reads_to_paths.csv missing from ${WINDOW}"
    fi
}


# Processing every shared window
while read WINDOW; do
    collectWindow ${WITH_SNPS_DIR}/output/run/${WINDOW} ${WITH_SNPS_DIR}/${LOG_CSV} ${WITH_SNPS_DIR}/${LOG_OPTIMIZER} ${WITH_SNPS_DIR}/${LOG_NPATHS} ${WITH_SNPS_DIR}/${LOG_NEDGES}
    collectWindow ${WITHOUT_SNPS_DIR}/output/run/${WINDOW} ${WITHOUT_SNPS_DIR}/${LOG_CSV} ${WITHOUT_SNPS_DIR}/${LOG_OPTIMIZER} ${WITHOUT_SNPS_DIR}/${LOG_NPATHS} ${WITHOUT_SNPS_DIR}/${LOG_NEDGES}    
done < shared_windows.txt
java OptimizerTimeToSeconds ${WITH_SNPS_DIR}/${LOG_OPTIMIZER} > ${WITH_SNPS_DIR}/log_optimizer_cleaned.txt
java OptimizerTimeToSeconds ${WITHOUT_SNPS_DIR}/${LOG_OPTIMIZER} > ${WITHOUT_SNPS_DIR}/log_optimizer_cleaned.txt
