#!/bin/bash
#

SAMPLE_IDS="1000151 1000513 1000920 1001399 1001980 1002322 1002826 1004266 1005038 1005444 1005938 1007198 1008775 1010384 1012440 1012736 1013536 1014457 1014625 1014694"
HPRC_MERGED_VCF="merged.regenotyped_kanpig_tp.dummy.split.vcf.gz"
TIME_COMMAND="time"
REFERENCE_FA="chm13v2.0.ebv.fa"
TRUVARI_BENCH_FLAGS="--sizemin 50 --sizemax 1000000 --sizefilt 50 --pctsize 0.9 --pctseq 0.9"
SV_LENGTHS="100 500 2500 10000 50000 100000"

set -euxo pipefail
LOG_FILE="results.log"
rm -rf ${LOG_FILE}
for ID in ${SAMPLE_IDS}; do
    rm -rf ./truvari/
    truvari bench ${TRUVARI_BENCH_FLAGS} --comp ${ID}.truvari_collapsed.vcf.gz --base ${HPRC_MERGED_VCF} --output ./truvari/
    TP=$(grep "\"TP-comp\":" ./truvari/summary.json | awk 'BEGIN {ORS=""} {print $2}')
    FP=$(grep "\"FP\":" ./truvari/summary.json | awk 'BEGIN {ORS=""} {print $2}')
    FN=$(grep "\"FN\":" ./truvari/summary.json | awk 'BEGIN {ORS=""} {print $2}')
    PR=$(grep "\"precision\":" ./truvari/summary.json | awk 'BEGIN {ORS=""} {print $2}')
    RE=$(grep "\"recall\":" ./truvari/summary.json | awk 'BEGIN {ORS=""} {print $2}')
    F1=$(grep "\"f1\":" ./truvari/summary.json | awk 'BEGIN {ORS=""} {print $2}')
    echo "${TP}${FP}${FN}${PR}${RE}${F1}" >> ${LOG_FILE}
    cat ${LOG_FILE}
    
    # Stratifying by SVLEN
    PREVIOUS_LENGTH="0"
    for SVLENGTH in ${SV_LENGTHS}; do
        FILTER_STRING="((SVLEN>${PREVIOUS_LENGTH} && SVLEN<=${SVLENGTH}) || (SVLEN>=-${SVLENGTH} && SVLEN<-${PREVIOUS_LENGTH}))"
        bcftools filter --include "${FILTER_STRING}" --output-type v ./truvari/tp-comp.vcf.gz > tmp1.vcf
        TP=$(bcftools view --no-header tmp1.vcf | wc -l)
        bcftools filter --include "${FILTER_STRING}" --output-type v ${ID}.truvari_collapsed.vcf.gz > tmp1.vcf
        TOTAL=$(bcftools view --no-header tmp1.vcf | wc -l)
        echo "${SVLENGTH},${TP},${TOTAL}" >> ${LOG_FILE}
        PREVIOUS_LENGTH=${SVLENGTH}
    done
    cat ${LOG_FILE}
done
