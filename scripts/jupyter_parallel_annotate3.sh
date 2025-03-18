#!/bin/bash
#
# Given VCFs that have already been annotated with TRs, annotates them using
# another BED.
#
SUBMISSION_ID=$1
NEW_BED_ADDRESS="??????????"
NEW_BED="???????????"
NEW_BED_GZ="${NEW_BED}.gz"
NEW_HDR_TXT="new.hdr.txt"

echo '##INFO=<ID=NEW,Number=1,Type=Integer,Description="90% of variant interval overlaps with the new BED">' > ${NEW_HDR_TXT}
gsutil cp ${NEW_BED_ADDRESS} .
awk -v OFS='\t' '{{ print $$0, 1 }}' ${NEW_BED} | /home/jupyter/callset_integration/edit/htslib-1.19.1/bgzip -c > ${NEW_BED_GZ}
/home/jupyter/callset_integration/edit/htslib-1.19.1/tabix -0 -s1 -b2 -e3 ${NEW_BED_GZ}

for VCF in $(ls scored-vcfs/${SUBMISSION_ID}/*score.tr.vcf.gz); do
    NEW_VCF=${VCF/score.tr.vcf.gz/score.tr.new.vcf.gz}
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools annotate --no-version -a ${NEW_BED_GZ} -h ${NEW_HDR_TXT} \
        -m +NEW_ONE -c CHROM,FROM,TO,- ${VCF} --min-overlap :0.9 | \
        grep -v '##INFO=<ID=NEW_ONE' | \
        sed -E 's/NEW_ONE/NEW=1/g' | \
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools annotate --no-version -a ${NEW_BED_GZ} -h ${NEW_HDR_TXT} \
        -m -NEW_ZERO -c CHROM,FROM,TO,- --min-overlap :0.9 | \
        grep -v '##INFO=<ID=NEW_ZERO' | \
        sed -E 's/NEW_ZERO/NEW=0/g' | \
    /home/jupyter/callset_integration/edit/bcftools-1.19/bcftools view -Oz -o ${NEW_VCF} &
done
wait
