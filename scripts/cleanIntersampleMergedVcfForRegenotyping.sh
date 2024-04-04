#!/bin/bash
#
# - Makes the inter-sample VCF a single-sample VCF with GT=0/1.
# - Replaces every call ID with a distinct integer.
#
ID=$1

set -euxo pipefail

rm -f ${ID}_cleaned.vcf*
./bcftools-1.18/bcftools view --header-only ${ID}.vcf.gz > ${ID}.txt
N_ROWS=$(wc -l < ${ID}.txt)
head -n $(( ${N_ROWS} - 1 )) ${ID}.txt > ${ID}_cleaned.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> ${ID}_cleaned.vcf
./bcftools-1.18/bcftools view --no-header ${ID}.vcf.gz | awk 'BEGIN{FS=OFS="\t"};{print $1, $2, i++, $4, $5, $6, $7, $8,"GT","0/1"}' >> ${ID}_cleaned.vcf
bgzip ${ID}_cleaned.vcf
tabix -f ${ID}_cleaned.vcf.gz