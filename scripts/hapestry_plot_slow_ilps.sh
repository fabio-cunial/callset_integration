#!/bin/bash
#
INPUT_DIR=$1   # Contains /output/run/<windows>
OUTPUT_DIR=$2
MAX_MINUTES="1"
WEIGHT_QUANTUM="1"
N_THREADS="1"
JAVA_FLAGS="-Xmx12G"

set -euxo pipefail


function collectWindow() {
    local WINDOW=$1
    local THREAD_ID=$2
    local STATS_FILE=${OUTPUT_DIR}/${THREAD_ID}.tsv
    local OUTPUT_DIR=${OUTPUT_DIR}/$(basename ${WINDOW})
    
    rm -rf ${OUTPUT_DIR}; mkdir ${OUTPUT_DIR}
    
    LOG_FILE="${WINDOW}/log.csv"
    if [ ! -e ${LOG_FILE} ]; then
        echo "WARNING: log.csv missing from ${WINDOW}"
        return 0
    fi
    
    grep 'feasibility,' ${LOG_FILE} > tmp_${THREAD_ID}.txt || echo ""
    N_LINES=$(wc -l < tmp_${THREAD_ID}.txt)
    if [ ${N_LINES} -eq 0 ]; then
        return 0
    fi
    HOURS_F=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 2 )
    MINUTES_F=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 3 )
    SUCCESS_F=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 6 )
    if [ ${SUCCESS_F} -eq 0 -o ${HOURS_F} -ge 1 -o ${MINUTES_F} -ge ${MAX_MINUTES} ]; then
        echo -e -n "1\t${HOURS_F}\t${MINUTES_F}\t${SUCCESS_F}\t" >> ${STATS_FILE}
        java ${JAVA_FLAGS} AnalyzeILPGraph ${WINDOW}/reads_to_paths.csv ${OUTPUT_DIR} ${WEIGHT_QUANTUM} 2>> ${STATS_FILE}
        return 0
    fi
    
    grep 'optimize_d,' ${LOG_FILE} > tmp_${THREAD_ID}.txt || echo ""
    N_LINES=$(wc -l < tmp_${THREAD_ID}.txt)
    if [ ${N_LINES} -eq 0 ]; then
        return 0
    fi
    HOURS_D=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 2 )
    MINUTES_D=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 3 )
    SUCCESS_D=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 6 )
    if [ ${SUCCESS_D} -eq 0 -o ${HOURS_D} -ge 1 -o ${MINUTES_D} -ge ${MAX_MINUTES} ]; then
        echo -e -n "2\t${HOURS_D}\t${MINUTES_D}\t${SUCCESS_D}\t" >> ${STATS_FILE}
        java ${JAVA_FLAGS} AnalyzeILPGraph ${WINDOW}/reads_to_paths.csv ${OUTPUT_DIR} ${WEIGHT_QUANTUM} 2>> ${STATS_FILE}
        return 0
    fi
    
    grep 'optimize_n_given_d,' ${LOG_FILE} > tmp_${THREAD_ID}.txt || echo ""
    N_LINES=$(wc -l < tmp_${THREAD_ID}.txt)
    if [ ${N_LINES} -eq 0 ]; then
        return 0
    fi
    HOURS_ND=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 2 )
    MINUTES_ND=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 3 )
    SUCCESS_ND=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 6 )
    if [ ${SUCCESS_ND} -eq 0 -o ${HOURS_ND} -ge 1 -o ${MINUTES_ND} -ge ${MAX_MINUTES} ]; then
        echo -e -n "3\t${HOURS_ND}\t${MINUTES_ND}\t${SUCCESS_ND}\t" >> ${STATS_FILE}
        java ${JAVA_FLAGS} AnalyzeILPGraph ${WINDOW}/reads_to_paths.csv ${OUTPUT_DIR} ${WEIGHT_QUANTUM} 2>> ${STATS_FILE}
        return 0
    fi
    
    grep 'optimize_d_plus_n,' ${LOG_FILE} > tmp_${THREAD_ID}.txt || echo ""
    N_LINES=$(wc -l < tmp_${THREAD_ID}.txt)
    if [ ${N_LINES} -eq 0 ]; then
        return
    fi
    HOURS_DN=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 2 )
    MINUTES_DN=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 3 )
    SUCCESS_DN=$( cat tmp_${THREAD_ID}.txt | cut -d , -f 6 )
    if [ ${SUCCESS_DN} -eq 0 -o ${HOURS_DN} -ge 1 -o ${MINUTES_DN} -ge ${MAX_MINUTES} ]; then
        echo -e -n "4\t${HOURS_DN}\t${MINUTES_DN}\t${SUCCESS_DN}\t" >> ${STATS_FILE}
        java ${JAVA_FLAGS} AnalyzeILPGraph ${WINDOW}/reads_to_paths.csv ${OUTPUT_DIR} ${WEIGHT_QUANTUM} 2>> ${STATS_FILE}
        return 0
    fi
}


function runThread() {
    local LIST=$1
    local THREAD_ID=$2
    
    while read WINDOW; do
        collectWindow ${INPUT_DIR}/output/run/${WINDOW} ${THREAD_ID}
    done < ${LIST}
}


# Main program
rm -rf ${OUTPUT_DIR}; mkdir ${OUTPUT_DIR}
ls ${INPUT_DIR}/output/run | sort > ${OUTPUT_DIR}/windows.txt
N_WINDOWS=$(wc -l < ${OUTPUT_DIR}/windows.txt)
echo "TOTAL WINDOWS: ${N_WINDOWS}"
split -l $(( ${N_WINDOWS} / ${N_THREADS} )) ${OUTPUT_DIR}/windows.txt prefix
i="0"
for LIST in $(ls prefix*); do
    runThread ${LIST} ${i} &
    i=$(( ${i} + 1 ))
done
wait
