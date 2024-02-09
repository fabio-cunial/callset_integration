#!/bin/bash
#

while read ADDRESS; do
    ID=$(basename ${ADDRESS} .vcf.gz)
    gsutil cp ${ADDRESS} .
    gunzip ${ID}.vcf.gz
    java PAV2SVs ${ID}.vcf 50 ${ID}_sv.vcf ${ID}_rest.vcf
    rm -f ${ID}_rest.vcf
    bgzip ${ID}_sv.vcf
    tabix ${ID}_sv.vcf.gz
done < list.txt