#!/bin/bash
INPUT_DIR=$1
TOOL=$2
CHROMOSOME=$3

set -euo pipefail

find ${INPUT_DIR}/chr${CHROMOSOME}_evaluation -type d -mindepth 1 -maxdepth 1 -print > directories_${CHROMOSOME}.txt
echo "Processing ${CHROMOSOME}"
while read DIRECTORY; do
    # Supported/unsupported VCFs
    if [ -e ${DIRECTORY}/${TOOL}/supported.vcf ]; then
        :
    else 
        echo "No supported.vcf in ${DIRECTORY}"
    fi
    if [ -e ${DIRECTORY}/${TOOL}/unsupported.vcf ]; then
        :
    else 
        echo "No unsupported.vcf in ${DIRECTORY}"
    fi
    
    # Successful execution
    SUCCESS=$(grep ${TOOL} ${DIRECTORY}/log.csv | cut -d , -f 6)
    if [ ${SUCCESS} -eq 0 ]; then
        echo "${TOOL} failed in ${DIRECTORY}"
    fi
done < directories_${CHROMOSOME}.txt
rm -f directories_${CHROMOSOME}.txt
