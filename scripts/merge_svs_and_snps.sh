#!/bin/bash
#
# Used for merging HPRC SNP VCFs (extracted from a joint SNP VCF) to existing
# SV VCFs.
#
LIST_FILE="list.txt"
SNP_VCF_DIR=".."
OUTPUT_DIR="svs_and_snps"

set -euxo pipefail


function merge() {
    local SAMPLE=$1
    local URL=$2
    
    gsutil cp ${URL} ${SAMPLE}.vcf.gz
    gsutil cp ${URL}.tbi ${SAMPLE}.vcf.gz.tbi
    bcftools concat --allow-overlaps --remove-duplicates --output-type z ${SAMPLE}.vcf.gz ${SNP_VCF_DIR}/${SAMPLE}.vcf.gz > ${OUTPUT_DIR}/${SAMPLE}_svs_and_snps.vcf.gz
    tabix -f ${OUTPUT_DIR}/${SAMPLE}_svs_and_snps.vcf.gz
}


rm -rf ${OUTPUT_DIR}; mkdir ${OUTPUT_DIR}
while read LINE; do
    SAMPLE=$(echo ${LINE} | cut -d , -f 1)
    URL=$(echo ${LINE} | cut -d , -f 2)
    merge ${SAMPLE} ${URL} &
done < ${LIST_FILE}
wait
