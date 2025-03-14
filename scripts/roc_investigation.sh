#!/bin/bash
#
TRUVARI_INTRASAMPLE_VCFS=$1
MIN_SV_LENGTH=$2

set -euxo pipefail

# Counts of raw calls per sample
PREFIX="${TRUVARI_INTRASAMPLE_VCFS}.${MIN_SV_LENGTH}"
N_INS="${PREFIX}.n_ins.txt"; rm -f ${N_INS}
N_DEL="${PREFIX}.n_del.txt"; rm -f ${N_DEL}
while read ADDRESS; do
    gsutil -m cp ${ADDRESS} ${PREFIX}.tmp1.vcf.gz
    gsutil -m cp ${ADDRESS}.tbi ${PREFIX}.tmp1.vcf.gz.tbi
    bcftools filter --include "SVTYPE=\"INS\" && (SVLEN<=-${MIN_SV_LENGTH} || SVLEN>=${MIN_SV_LENGTH})" --output-type z ${PREFIX}.tmp1.vcf.gz > ${PREFIX}.tmp2.vcf.gz
    tabix -f ${PREFIX}.tmp2.vcf.gz
    bcftools view --no-header ${PREFIX}.tmp2.vcf.gz | wc -l >> ${N_INS}
    bcftools filter --include "SVTYPE=\"DEL\" && (SVLEN<=-${MIN_SV_LENGTH} || SVLEN>=${MIN_SV_LENGTH})" --output-type z ${PREFIX}.tmp1.vcf.gz > ${PREFIX}.tmp2.vcf.gz
    tabix -f ${PREFIX}.tmp2.vcf.gz
    bcftools view --no-header ${PREFIX}.tmp2.vcf.gz | wc -l >> ${N_DEL}
done < ${TRUVARI_INTRASAMPLE_VCFS}
