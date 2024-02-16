#!/bin/bash
#
MIN_SV_LENGTH="20";

set -euxo pipefail


while read LINE; do
    ID=${LINE%,*}
    ADDRESS=${LINE#*,}
    gsutil cp ${ADDRESS} ./${ID}.vcf.gz
    tabix ${ID}.vcf.gz
    
    # Making sure that the VCF has the right sample name
    echo ${ID} > samples.txt
    bcftools reheader --samples samples.txt ${ID}.vcf.gz > reheaded.vcf.gz
    tabix reheaded.vcf.gz
    rm -f ${ID}.vcf.gz*

    # Making sure that the VCF has no multiallelic records
    bcftools norm --multiallelics - --output-type v reheaded.vcf.gz > fixed.vcf
    rm -f reheaded.vcf.gz*
    
    java Dipcall2VCF fixed.vcf ${MIN_SV_LENGTH} ${ID}_sv.vcf
    bgzip ${ID}_sv.vcf
    tabix ${ID}_sv.vcf.gz
done < list.txt
