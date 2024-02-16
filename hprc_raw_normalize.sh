#!/bin/bash
#
ROOT_DIR="gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HPRC"

set -euxo pipefail

while read ID; do
    gsutil cp ${ROOT_DIR}/${ID}/truth_chm13_raw_sv.vcf.gz .
    bcftools norm -m - truth_chm13_raw_sv.vcf.gz > tmp.vcf
    java HPRCFilter tmp.vcf . 1
    rm -f tmp.vcf
    bgzip ${ID}.vcf
    tabix ${ID}.vcf.gz
    gsutil mv ${ID}.vcf.gz ${ROOT_DIR}/${ID}/truth_chm13_raw_sv_norm.vcf.gz
    gsutil mv ${ID}.vcf.gz.tbi ${ROOT_DIR}/${ID}/truth_chm13_raw_sv_norm.vcf.gz.tbi
done < hprc_samples.txt
