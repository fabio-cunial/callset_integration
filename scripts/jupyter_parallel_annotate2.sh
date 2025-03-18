#!/bin/bash
#
# Subsets to a confident BED and then annotates TRs.
#
SUBMISSION_ID=$1
TR_BED="human_GRCh38_no_alt_analysis_set.trf.bed"
TR_BED_GZ="${TR_BED}.gz"
TR_HDR_TXT="tr.hdr.txt"
CONFIDENT_BED="GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed"

echo '##INFO=<ID=TR,Number=1,Type=Integer,Description="90% of variant interval overlaps with simple-repeat region">' > ${TR_HDR_TXT}
gsutil cp gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/references/${TR_BED} .
awk -v OFS='\t' '{{ print $$0, 1 }}' ${TR_BED} | /home/jupyter/callset_integration/edit/htslib-1.19.1/bgzip -c > ${TR_BED_GZ}
/home/jupyter/callset_integration/edit/htslib-1.19.1/tabix -0 -s1 -b2 -e3 ${TR_BED_GZ}


function annotate() {
    local VCF=$1
    local i=$2
    local TR_VCF=${VCF/score.vcf.gz/score.tr.vcf.gz}
    
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools view --regions-file ${CONFIDENT_BED} --regions-overlap pos --output-type z ${VCF} > tmp_${i}.vcf.gz
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools index --force tmp_${i}.vcf.gz
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools annotate --no-version -a ${TR_BED_GZ} -h ${TR_HDR_TXT} \
        -m +TR_ONE -c CHROM,FROM,TO,- tmp_${i}.vcf.gz --min-overlap :0.9 tmp_${i}.vcf.gz | \
        grep -v '##INFO=<ID=TR_ONE' | \
        sed -E 's/TR_ONE/TR=1/g' | \
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools annotate --no-version -a ${TR_BED_GZ} -h ${TR_HDR_TXT} \
        -m -TR_ZERO -c CHROM,FROM,TO,- --min-overlap :0.9 | \
        grep -v '##INFO=<ID=TR_ZERO' | \
        sed -E 's/TR_ZERO/TR=0/g' | \
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools view -Oz -o ${TR_VCF}
    rm -f tmp_${i}.vcf.gz*
}


i="0"
for VCF in $(ls scored-vcfs/${SUBMISSION_ID}/*score.vcf.gz); do
    annotate ${VCF} ${i} &
    i=$(( ${i} + 1 ))
done
wait
