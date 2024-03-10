#!/bin/bash
#
INPUT_DIR="gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/HGSVC2"
OUTPUT_DIR="gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/bcftools_merge_dipcall/hgsvc2/new_per_sample_dipcall_vcfs"
SAMPLES="HG00512 HG00513 HG00514 HG00731 HG00732 HG00733 HG02818 HG03125 HG03486 NA12878 NA19238 NA19239 NA19240 NA24385"
MIN_SV_LENGTH="20"

set -euxo pipefail

for SAMPLE in ${SAMPLES}; do
    gsutil -m cp ${INPUT_DIR}/${SAMPLE}/${SAMPLE}.dipcall.vcf.'gz*' .
    bcftools norm --multiallelics - --output-type v ${SAMPLE}.dipcall.vcf.gz > tmp.vcf
    rm -f ${SAMPLE}.dipcall.vcf.gz
    java Dipcall2VCF tmp.vcf ${MIN_SV_LENGTH} ${SAMPLE}.dipcall_sv.vcf
    rm -f tmp.vcf
    bgzip ${SAMPLE}.dipcall_sv.vcf
    tabix ${SAMPLE}.dipcall_sv.vcf.gz
    gsutil -m mv ${SAMPLE}.dipcall_sv.vcf.'gz*' ${OUTPUT_DIR}
done
