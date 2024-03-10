#!/bin/bash
INPUT_DIR=$1

set -euo pipefail

TOOL="sniffles"
SUPPORTED_VCF_FILE="${INPUT_DIR}/supported_all.vcf"; rm -f ${SUPPORTED_VCF_FILE}
UNSUPPORTED_VCF_FILE="${INPUT_DIR}/unsupported_all.vcf"; rm -f ${UNSUPPORTED_VCF_FILE}

# Collecting
for CHROMOSOME in $(seq 11 22); do
    find ${INPUT_DIR}/chr${CHROMOSOME}_evaluation -type d -mindepth 1 -maxdepth 1 -print > directories.txt
    while read DIRECTORY; do
        if [ -e ${DIRECTORY}/${TOOL}/supported.vcf ]; then
            SUPPORTED=$(wc -l < ${DIRECTORY}/${TOOL}/supported.vcf)
            SUPPORTED=$(( ${SUPPORTED} - 1 ))
            tail -n +2 ${DIRECTORY}/${TOOL}/supported.vcf >> ${SUPPORTED_VCF_FILE}
        else
            echo "No supported file in ${DIRECTORY}"
        fi
        if [ -e ${DIRECTORY}/${TOOL}/unsupported.vcf ]; then
            UNSUPPORTED=$(wc -l < ${DIRECTORY}/${TOOL}/unsupported.vcf)
            UNSUPPORTED=$(( ${UNSUPPORTED} - 1 ))
            tail -n +2 ${DIRECTORY}/${TOOL}/unsupported.vcf >> ${UNSUPPORTED_VCF_FILE}
        else
            echo "No unsupported file in ${DIRECTORY}"
        fi
    done < directories.txt
done

# Counting
SUPPORTED_DEL=$(grep 'SVTYPE=DEL' ${SUPPORTED_VCF_FILE} | wc -l)
SUPPORTED_INS=$(grep 'SVTYPE=INS' ${SUPPORTED_VCF_FILE} | wc -l)
UNSUPPORTED_DEL=$(grep 'SVTYPE=DEL' ${UNSUPPORTED_VCF_FILE} | wc -l)
UNSUPPORTED_INS=$(grep 'SVTYPE=INS' ${UNSUPPORTED_VCF_FILE} | wc -l)
echo "Supported DEL=${SUPPORTED_DEL} Unsupported DEL=${UNSUPPORTED_DEL}"   # 27%
echo "Supported INS=${SUPPORTED_INS} Unsupported INS=${UNSUPPORTED_INS}"   # 42%

# UNSUPPORTED_PBSV_FILE="${INPUT_DIR}/unsupported_all_pbsv.vcf"
# cut -f 3,8 ${UNSUPPORTED_VCF_FILE} | grep pbsv > ${UNSUPPORTED_PBSV_FILE}
# UNSUPPORTED_DEL=$(grep 'SVTYPE=DEL' ${UNSUPPORTED_PBSV_FILE} | wc -l)
# UNSUPPORTED_INS=$(grep 'SVTYPE=INS' ${UNSUPPORTED_PBSV_FILE} | wc -l)
# echo "Unsupported DEL from PBSV: ${UNSUPPORTED_DEL}   Unsupported INS from PBSV: ${UNSUPPORTED_INS}"
#
# UNSUPPORTED_SNIFFLES_FILE="${INPUT_DIR}/unsupported_all_sniffles.vcf"
# cut -f 3,8 ${UNSUPPORTED_VCF_FILE} | grep Sniffles > ${UNSUPPORTED_SNIFFLES_FILE}
# UNSUPPORTED_DEL=$(grep 'SVTYPE=DEL' ${UNSUPPORTED_SNIFFLES_FILE} | wc -l)
# UNSUPPORTED_INS=$(grep 'SVTYPE=INS' ${UNSUPPORTED_SNIFFLES_FILE} | wc -l)
# echo "Unsupported DEL from SNIFFLES: ${UNSUPPORTED_DEL}   Unsupported INS from SNIFFLES: ${UNSUPPORTED_INS}"
