#!/bin/bash
#
MIN_SV_LENGTH="10"
INPUT_DIR="gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC"
OUTPUT_DIR="gs://fc-42d72f6a-6be3-4ced-888e-06d4844681c3/dipcall_${MIN_SV_LENGTH}bp"
FULL_DIPCALL_FILE="/Users/fcunial/git/callset_integration/21309.txt"  # "/Users/fcunial/git/callset_integration/HPRC_chm13_full_dipcall_from_infogain.txt"
N_THREADS="8";

set -euxo pipefail

while read ROW; do
    SAMPLE=$(echo ${ROW} | cut -d ' ' -f 1); ADDRESS=$(echo ${ROW} | cut -d ' ' -f 2)
    gsutil cp ${ADDRESS} tmp1.vcf.gz
    bcftools norm --threads ${N_THREADS} --multiallelics - --output-type v tmp1.vcf.gz > tmp2.vcf
    rm -f tmp1.vcf.gz
    java Dipcall2VCF tmp2.vcf ${MIN_SV_LENGTH} ${SAMPLE}.dipcall_${MIN_SV_LENGTH}bp.vcf
    rm -f tmp2.vcf
    bgzip -@ ${N_THREADS} ${SAMPLE}.dipcall_${MIN_SV_LENGTH}bp.vcf
    tabix -f ${SAMPLE}.dipcall_${MIN_SV_LENGTH}bp.vcf.gz
    gsutil mv ${SAMPLE}.dipcall_${MIN_SV_LENGTH}bp.vcf.'gz*' ${OUTPUT_DIR}
done < ${FULL_DIPCALL_FILE}
