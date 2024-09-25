#!/bin/bash
#
# Collects all the records in the supported VCFs created by hapestry evaluate.
#
ROOT_DIR=$1
SUPPORTED_VCF=$2
TOOL_ID="SVs_and_SNPs"

rm -f tmp.vcf
ls -d ${ROOT_DIR}/* > list.txt
while read DIRECTORY; do
    FILE="${DIRECTORY}/${TOOL_ID}/supported.vcf"
    if [ -e ${FILE} ]; then
        tail -n +3 ${FILE} >> tmp.vcf
    fi
done < list.txt
rm -f list.txt
sort tmp.vcf > ${SUPPORTED_VCF}
rm -f tmp.vcf
