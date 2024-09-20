#!/bin/bash
#
# Used for creating SNP VCFs for each HPRC sample, given the joint-called SNP
# VCF of all HPRC samples.
#
# Remark: this implementation is naive and should be improved by using bcftools
# to split and select nonzero GTs (instead of grep).
#
SAMPLES="HG002 HG00438 HG005 HG00621 HG00673 HG00733 HG00735 HG00741 HG01071 HG01106 HG01109 HG01123 HG01175 HG01243 HG01258  HG01358 HG01361 HG01891 HG01928 HG01952 HG01978 HG02055 HG02080 HG02109 HG02145 HG02148 HG02257 HG02486 HG02559 HG02572 HG02622 HG02630 HG02717 HG02723 HG02818 HG02886 HG03098 HG03453 HG03486 HG03492 HG03516 HG03540 HG03579 NA18906 NA19240 NA20129 NA21309"
JOINT_VCF="HPRC.DV.joint.g.vcf.bgz"

set -euxo pipefail


function project() {
    local SAMPLE=$1
    local HEADER_FILE=$2
    
    rm -f ${SAMPLE}.vcf*
    cp ${HEADER_FILE} ${SAMPLE}.vcf
    echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${SAMPLE}" >> ${SAMPLE}.vcf
    bcftools view --no-header --samples ${SAMPLE} --output-type v ${JOINT_VCF} | grep '0/1\|1/0\|1/1' >> ${SAMPLE}.vcf
    bgzip ${SAMPLE}.vcf
    tabix -f ${SAMPLE}.vcf.gz
}


bcftools view --header-only ${JOINT_VCF} > tmp.txt
N_ROWS=$(wc -l < tmp.txt)
head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
rm -f tmp.txt
for SAMPLE in ${SAMPLES}; do
    project ${SAMPLE} header.txt &
done
wait
