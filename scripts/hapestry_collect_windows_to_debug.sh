#!/bin/bash
#
#
WINDOWS_FILE="windows.txt"

set -euxo pipefail

# Very slow, can be exectued just once.
rm -f tmp.txt
find ./hapestry_merge/ -type d -exec basename {} ';' > tmp.txt
tail -n +4 tmp.txt | tr '_' '\t' | tr '-' '\t' > hapestry_merge_windows.txt
rm -f tmp.txt

while read WINDOW; do
    # Evaluation window
    rm -f ${WINDOW}_evaluation.tar.gz
    tar -czf ${WINDOW}_evaluation.tar.gz chr1_evaluation/${WINDOW}
    
    # Contained hapestry merge windows
    TOKENS=$(echo ${WINDOW} | tr '_' ',' | tr '-' ',')
    WINDOW_START=$(echo ${TOKENS} | cut -d , -f 2)
    WINDOW_END=$(echo ${TOKENS} | cut -d , -f 3)
    awk -v start=${WINDOW_START} -v end=${WINDOW_END} '{ if ($1=="chr1" && $2>=start && $3<=end) print $1 "," $2 "," $3 }' hapestry_merge_windows.txt > contained.txt
    if [ $(wc -l < contained.txt) -gt 0 ]; then
        while read WINDOW2; do
            WINDOW_START=$(echo ${WINDOW2} | cut -d , -f 2)
            WINDOW_END=$(echo ${WINDOW2} | cut -d , -f 3)
            ARCHIVE_NAME="${WINDOW}_hapestrymerge_${WINDOW_START}_${WINDOW_END}.tar.gz"
            rm -f ${ARCHIVE_NAME}
            tar -czf ${ARCHIVE_NAME} $(find ./hapestry_merge/ -type d -name chr1_${WINDOW_START}-${WINDOW_END})
        done < contained.txt
    else
        # Containing hapestry merge windows
        awk -v start=${WINDOW_START} -v end=${WINDOW_END} '{ if ($1=="chr1" && $2<=start && $3>=end) print $1 "," $2 "," $3 }' hapestry_merge_windows.txt > containing.txt
        #printf("%s,%d,%d",$1,$2,$3)
        if [ $(wc -l < containing.txt) -gt 0 ]; then
            while read WINDOW2; do
                WINDOW_START=$(echo ${WINDOW2} | cut -d , -f 2)
                WINDOW_END=$(echo ${WINDOW2} | cut -d , -f 3)
                ARCHIVE_NAME="${WINDOW}_hapestrymerge_${WINDOW_START}_${WINDOW_END}.tar.gz"
                rm -f ${ARCHIVE_NAME}
                tar -czf ${ARCHIVE_NAME} $(find ./hapestry_merge/ -type d -name chr1_${WINDOW_START}-${WINDOW_END})
            done < containing.txt
        else
            echo "ERROR: No hapestry merge window is contained in or contains window ${WINDOW}"
        fi
    fi
done < ${WINDOWS_FILE}
