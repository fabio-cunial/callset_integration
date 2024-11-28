#!/bin/bash
#
# Remark: the script assumes an intra-sample-merged VCF in input, containing
# calls from pbsv, sniffles, pav.
#
# SAMPLE_ID="q100"
ROOT_DIR="/Users/fcunial/Downloads/sam_ks"
CALLER_VCF="${ROOT_DIR}/original_regenotyped.vcf.gz"
# ALIGNMENTS_BAM="${ROOT_DIR}/HG002.bam"
# READS_FASTQ_GZ="${ROOT_DIR}/HG002.fastq.gz"
TRUTH_VCF="${ROOT_DIR}/original_truth.vcf.gz"
TRUTH_BED="${ROOT_DIR}/original_truth.bed"
REFERENCE_FA="/Users/fcunial/Downloads/sniffles-gt/GCA_000001405.15_GRCh38_no_alt_analysis_set.fa"
TR_BED="/Users/fcunial/Downloads/sniffles-gt/GRCh38_AllTandemRepeatsandHomopolymers_slop5.bed"
NONTR_BED="/Users/fcunial/Downloads/sniffles-gt/GRCh38_notinAllTandemRepeatsandHomopolymers_slop5.bed"

TR_SMALL_OVERLAP="0.1"
TR_LARGE_OVERLAP="0.9"
SVLEN_BINS="50,100,500,1000,5000,10000,50000";
SVLEN_MIN="0"; SVLEN_MAX="50000";
TRUVARI_BENCH_SETTINGS="--sizefilt ${SVLEN_MIN} --sizemax ${SVLEN_MAX} --sizemin 0"

RESOLVE_COMMAND="python /Users/fcunial/git/long-reads-public-codebase/docker/truvari_intrasample/resolve.py"
# KANPIG_COMMAND="/Users/fcunial/Downloads/kanpig"
CLEAN_VCF_COMMAND="java CleanVCF"
N_THREADS="8";


# Makes sure that the truth and merged VCF are in the right format and contain
# only a specific set of calls.
#
# @param 4: keep only calls with: 0=no filter; 1=large overlap with the TR
# track; 2=small overlap with the TR track; 3=large overlap with the non-TR
# track.
#
function formatVcf() {
    local INPUT_VCF_GZ=$1
    local OUTPUT_VCF_GZ=$2
    local IS_TRUTH=$3
    local TR_STATUS=$4
    
    # Removing multiallelic records
    rm -f tmp0.vcf.gz*
    bcftools norm --multiallelics - --output-type z ${INPUT_VCF_GZ} > tmp0.vcf.gz
    tabix -f tmp0.vcf.gz

    # Adding REF and ALT with Adam's script, and removing BNDs.
    rm -f ${ROOT_DIR}/tmp1.vcf*
    ${RESOLVE_COMMAND} tmp0.vcf.gz ${REFERENCE_FA} | bcftools view -i "SVTYPE!='BND'" | bcftools sort -O z -o ${ROOT_DIR}/tmp1.vcf.gz
    tabix -f ${ROOT_DIR}/tmp1.vcf.gz

    # Keeping only calls in the truth BED (any overlap).
    rm -rf ${ROOT_DIR}/tmp2.vcf*
    bcftools view --header-only ${ROOT_DIR}/tmp1.vcf.gz > ${ROOT_DIR}/tmp2.vcf
    bedtools intersect -u -a ${ROOT_DIR}/tmp1.vcf.gz -b ${TRUTH_BED} >> ${ROOT_DIR}/tmp2.vcf
    rm -f ${ROOT_DIR}/tmp2.vcf.gz*; bgzip ${ROOT_DIR}/tmp2.vcf
    tabix -f ${ROOT_DIR}/tmp2.vcf.gz

    # If not a truth VCF: keeping only calls in the given TR track (with min.
    # overlap).
    if [ ${IS_TRUTH} -eq 1 ]; then
        cp ${ROOT_DIR}/tmp2.vcf.gz ${ROOT_DIR}/tmp3.vcf.gz
        cp ${ROOT_DIR}/tmp2.vcf.gz.tbi ${ROOT_DIR}/tmp3.vcf.gz.tbi
    else
        rm -rf ${ROOT_DIR}/tmp3.vcf*
        if [ ${TR_STATUS} -eq 0 ]; then
            cp ${ROOT_DIR}/tmp2.vcf.gz ${ROOT_DIR}/tmp3.vcf.gz
            cp ${ROOT_DIR}/tmp2.vcf.gz.tbi ${ROOT_DIR}/tmp3.vcf.gz.tbi
        elif [ ${TR_STATUS} -eq 1 ]; then
            bcftools view --header-only ${ROOT_DIR}/tmp2.vcf.gz > ${ROOT_DIR}/tmp3.vcf
            bedtools intersect -u -f ${TR_LARGE_OVERLAP} -a ${ROOT_DIR}/tmp2.vcf.gz -b ${TR_BED} >> ${ROOT_DIR}/tmp3.vcf
            rm -f ${ROOT_DIR}/tmp3.vcf.gz*; bgzip ${ROOT_DIR}/tmp3.vcf
            tabix -f ${ROOT_DIR}/tmp3.vcf.gz
        elif [ ${TR_STATUS} -eq 2 ]; then
            bcftools view --header-only ${ROOT_DIR}/tmp2.vcf.gz > ${ROOT_DIR}/tmp3.vcf
            bedtools intersect -u -f ${TR_SMALL_OVERLAP} -a ${ROOT_DIR}/tmp2.vcf.gz -b ${TR_BED} >> ${ROOT_DIR}/tmp3.vcf
            rm -f ${ROOT_DIR}/tmp3.vcf.gz*; bgzip ${ROOT_DIR}/tmp3.vcf
            tabix -f ${ROOT_DIR}/tmp3.vcf.gz
        elif [ ${TR_STATUS} -eq 3 ]; then
            bcftools view --header-only ${ROOT_DIR}/tmp2.vcf.gz > ${ROOT_DIR}/tmp3.vcf
            bedtools intersect -u -f ${TR_LARGE_OVERLAP} -a ${ROOT_DIR}/tmp2.vcf.gz -b ${NONTR_BED} >> ${ROOT_DIR}/tmp3.vcf
            rm -f ${ROOT_DIR}/tmp3.vcf.gz*; bgzip ${ROOT_DIR}/tmp3.vcf
            tabix -f ${ROOT_DIR}/tmp3.vcf.gz
        else
            :
        fi
    fi

    # If not a truth VCF: keeping only calls in the given length range.
    if [ ${IS_TRUTH} -eq 1 ]; then
        cp ${ROOT_DIR}/tmp3.vcf.gz ${ROOT_DIR}/tmp4.vcf.gz
        cp ${ROOT_DIR}/tmp3.vcf.gz.tbi ${ROOT_DIR}/tmp4.vcf.gz.tbi
    else
        rm -rf ${ROOT_DIR}/tmp4.vcf*
        FILTER_STRING="SVLEN>=${SVLEN_MIN} && SVLEN<=${SVLEN_MAX}"
        bcftools filter -i "${FILTER_STRING}" --output-type z ${ROOT_DIR}/tmp3.vcf.gz > ${ROOT_DIR}/tmp4.vcf.gz
        tabix ${ROOT_DIR}/tmp4.vcf.gz
    fi

    # Finalizing
    cp ${ROOT_DIR}/tmp4.vcf.gz ${OUTPUT_VCF_GZ}
    cp ${ROOT_DIR}/tmp4.vcf.gz.tbi ${OUTPUT_VCF_GZ}.tbi
    rm -f ${ROOT_DIR}/tmp*.vcf*
}


# @param 3: 0=no filter; 1=large overlap with the TR track; 2=small overlap
# with the TR track; 3=large overlap with the non-TR track.
#
function afterRegenotyping() {
    local TR_STATUS=$1
    local GENOTYPER=$2
    local LOG_FILE=$3
    
    tabix -f ${ROOT_DIR}/regenotyped.vcf.gz
    N_CALLS_00_AFTER=$(bcftools filter -i "${FILTER_STRING}" ${ROOT_DIR}/regenotyped.vcf.gz | grep '^[^#]' | wc -l || echo 0)
    VALUE_BEFORE="${N_CALLS_00_BEFORE}/${N_CALLS}"; 
    #VALUE_BEFORE=$(bc -l <<< ${VALUE_BEFORE})
    VALUE_AFTER="${N_CALLS_00_AFTER}/${N_CALLS}"; 
    #VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
    echo "${TR_STATUS},${GENOTYPER} Fraction of 0/0 calls: Before regenotyping: ${VALUE_BEFORE} -> After regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
    if [ ${GENOTYPER} = sniffles -o ${GENOTYPER} = cutesv ]; then
        N_CALLS_DV_ZERO_AFTER=$(bcftools filter -i 'DV=0' ${ROOT_DIR}/regenotyped.vcf.gz | grep '^[^#]' | wc -l)
        VALUE_AFTER="${N_CALLS_DV_ZERO_AFTER}/${N_CALLS}"; 
        #VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
        echo "${TR_STATUS},${GENOTYPER} Fraction of calls with 0 supporting reads after regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
        bcftools view ${ROOT_DIR}/merged.vcf.gz > ${ROOT_DIR}/merged.vcf
        bcftools view ${ROOT_DIR}/regenotyped.vcf.gz > ${ROOT_DIR}/regenotyped.vcf
        java SupportedByZeroReads DV ${ROOT_DIR}/merged.vcf ${ROOT_DIR}/regenotyped.vcf ${SVLEN_MIN} ${SVLEN_MAX} ${SVLEN_BINS} > ${ROOT_DIR}/${TR_STATUS}_${GENOTYPER}_zeroReads.log
        rm -f ${ROOT_DIR}/merged.vcf ${ROOT_DIR}/regenotyped.vcf
    elif [ ${GENOTYPER} = kanpig ]; then
        # N_CALLS_DV_ZERO_AFTER=$(bcftools filter -i 'FORMAT/AD[0:1]=0' ${ROOT_DIR}/regenotyped.vcf.gz | grep '^[^#]' | wc -l)
        # VALUE_AFTER="${N_CALLS_DV_ZERO_AFTER}/${N_CALLS}";
        #VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
        # echo "${TR_STATUS},${GENOTYPER} Fraction of calls with 0 supporting reads after regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
        # bcftools view ${ROOT_DIR}/merged.vcf.gz > ${ROOT_DIR}/merged.vcf
        # bcftools view ${ROOT_DIR}/regenotyped.vcf.gz > ${ROOT_DIR}/regenotyped.vcf
        # java SupportedByZeroReads AD ${ROOT_DIR}/merged.vcf ${ROOT_DIR}/regenotyped.vcf ${SVLEN_MIN} ${SVLEN_MAX} ${SVLEN_BINS} > ${ROOT_DIR}/${TR_STATUS}_${GENOTYPER}_zeroReads.log
        # rm -f ${ROOT_DIR}/merged.vcf ${ROOT_DIR}/regenotyped.vcf
        :
    elif [ ${GENOTYPER} = svjedigraph ]; then
        N_CALLS_DV_ZERO_AFTER=$(bcftools filter -i 'FORMAT/AD[0:1]="0"' ${ROOT_DIR}/regenotyped.vcf.gz | grep '^[^#]' | wc -l)
        VALUE_AFTER="${N_CALLS_DV_ZERO_AFTER}/${N_CALLS}"; 
        #VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
        echo "${TR_STATUS},${GENOTYPER} Fraction of calls with 0 supporting reads after regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
        bcftools view ${ROOT_DIR}/merged.vcf.gz > ${ROOT_DIR}/merged.vcf
        bcftools view ${ROOT_DIR}/regenotyped.vcf.gz > ${ROOT_DIR}/regenotyped.vcf
        java SupportedByZeroReads AD ${ROOT_DIR}/merged.vcf ${ROOT_DIR}/regenotyped.vcf ${SVLEN_MIN} ${SVLEN_MAX} ${SVLEN_BINS} > ${ROOT_DIR}/${TR_STATUS}_${GENOTYPER}_zeroReads.log
        rm -f ${ROOT_DIR}/merged.vcf ${ROOT_DIR}/regenotyped.vcf
    elif [ ${GENOTYPER} = sams ]; then
        :
    fi

    # Computing TPs and ROC curves
    rm -rf ${ROOT_DIR}/truvari/
    truvari bench ${TRUVARI_BENCH_SETTINGS} -b ${ROOT_DIR}/truth.vcf.gz -c ${ROOT_DIR}/regenotyped.vcf.gz -o ${ROOT_DIR}/truvari/
    bcftools view ${ROOT_DIR}/regenotyped.vcf.gz > ${ROOT_DIR}/regenotyped.vcf
    bcftools view ${ROOT_DIR}/truvari/tp-comp.vcf.gz > ${ROOT_DIR}/truvari/tp-comp.vcf
    bcftools view ${ROOT_DIR}/truth.vcf.gz > ${ROOT_DIR}/truth.vcf
    if [ ${GENOTYPER} = sniffles -o ${GENOTYPER} = cutesv ]; then
        java ROCcurve GQ 1 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        java ROCcurve GQ 0 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        java ROCcurve DV 1 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        java ROCcurve DV 0 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        bcftools filter -i 'DV=0' ${ROOT_DIR}/regenotyped.vcf.gz > ${ROOT_DIR}/${TR_STATUS}_${GENOTYPER}_zeroReads.vcf
    elif [ ${GENOTYPER} = kanpig ]; then
        java ROCcurve KS 1 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        # java ROCcurve SQ 0 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        # java ROCcurve AD 1 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        # java ROCcurve AD 0 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        # bcftools filter -i 'FORMAT/AD[0:1]=0' ${ROOT_DIR}/regenotyped.vcf.gz > ${ROOT_DIR}/${TR_STATUS}_${GENOTYPER}_zeroReads.vcf
    elif [ ${GENOTYPER} = svjedigraph ]; then
        #java ROCcurve PL 1 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        #java ROCcurve PL 0 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        java ROCcurve AD 1 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        java ROCcurve AD 0 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        bcftools filter -i 'FORMAT/AD[0:1]="0"' ${ROOT_DIR}/regenotyped.vcf.gz > ${ROOT_DIR}/${TR_STATUS}_${GENOTYPER}_zeroReads.vcf
    elif [ ${GENOTYPER} = sams ]; then
        java ROCcurve SCORE 1 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
        java ROCcurve SCORE 0 ${SVLEN_MIN} ${SVLEN_MAX} ${ROOT_DIR}/regenotyped.vcf ${ROOT_DIR}/truvari/tp-comp.vcf ${ROOT_DIR}/truth.vcf ${ROOT_DIR} ${TR_STATUS} ${GENOTYPER} ${SVLEN_BINS}
    fi
}




# Main program
set -euxo pipefail

rm -f *.log
for TR_STATUS in 0 1 2 3; do
    LOG_FILE="${ROOT_DIR}/${TR_STATUS}.log"
    rm -f ${LOG_FILE}

    # Formatting truth and merged VCF
    formatVcf ${CALLER_VCF} ${ROOT_DIR}/merged.vcf.gz 0 ${TR_STATUS}
    truvari anno svinfo --minsize 0 ${TRUTH_VCF} | bgzip > ${ROOT_DIR}/tmp-1.vcf.gz
    tabix ${ROOT_DIR}/tmp-1.vcf.gz
    formatVcf ${ROOT_DIR}/tmp-1.vcf.gz ${ROOT_DIR}/truth.vcf.gz 1 ${TR_STATUS}

    # Regenotyping
    rm -f ${ROOT_DIR}/regenotyped.vcf.gz
    N_CALLS=$(bcftools view --no-header ${ROOT_DIR}/merged.vcf.gz | wc -l)
    FILTER_STRING="GT=\"./.\" || GT=\"./0\" || GT=\"0/.\" || GT=\"0/0\" || GT=\".|.\" || GT=\".|0\" || GT=\"0|.\" || GT=\"0|0\""
    N_CALLS_00_BEFORE=$(bcftools filter -i "${FILTER_STRING}" ${ROOT_DIR}/merged.vcf.gz | grep '^[^#]' | wc -l || echo 0)

    # # SNIFFLES FORCE
#     sniffles --threads ${N_THREADS} --reference ${REFERENCE_FA} --input ${ALIGNMENTS_BAM} --genotype-vcf ${ROOT_DIR}/merged.vcf.gz --vcf ${ROOT_DIR}/regenotyped.vcf.gz
#     afterRegenotyping ${TR_STATUS} sniffles ${LOG_FILE}

    # # CUTESV FORCE
#     rm -rf ./cutesv_tmp ; mkdir ./cutesv_tmp
#     cuteSV --threads ${N_THREADS} -Ivcf ${ROOT_DIR}/merged.vcf.gz --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 1 --max_size -1 ${ALIGNMENTS_BAM} ${REFERENCE_FA} ${ROOT_DIR}/tmp1.vcf.gz ./cutesv_tmp
#     bcftools sort --output-type z ${ROOT_DIR}/tmp1.vcf.gz > ${ROOT_DIR}/regenotyped.vcf.gz
#     rm -f ${ROOT_DIR}/tmp1.vcf.gz
#     afterRegenotyping ${TR_STATUS} cutesv ${LOG_FILE}

    # KANPIG
    # ${KANPIG_COMMAND} --threads ${N_THREADS} --sizemin ${SVLEN_MIN} --sizemax ${SVLEN_MAX} --input ${ROOT_DIR}/merged.vcf.gz --bam ${ALIGNMENTS_BAM} --reference ${REFERENCE_FA} --out ${ROOT_DIR}/tmp1.vcf.gz
    # bcftools sort --output-type z ${ROOT_DIR}/tmp1.vcf.gz > ${ROOT_DIR}/regenotyped.vcf.gz
    # rm -f ${ROOT_DIR}/tmp1.vcf.gz
    
    cp ${ROOT_DIR}/merged.vcf.gz ${ROOT_DIR}/regenotyped.vcf.gz
    cp ${ROOT_DIR}/merged.vcf.gz.tbi ${ROOT_DIR}/regenotyped.vcf.gz.tbi
    
    afterRegenotyping ${TR_STATUS} kanpig ${LOG_FILE}

    # LRCALLER
    #
    #

    # Callers supporting a call
    # rm -rf ${ROOT_DIR}/truvari/
#     truvari bench ${TRUVARI_BENCH_SETTINGS} -b ${ROOT_DIR}/truth.vcf.gz -c ${ROOT_DIR}/merged.vcf.gz -o ${ROOT_DIR}/truvari/
#     bcftools view ${ROOT_DIR}/truvari/tp-comp.vcf.gz > ${ROOT_DIR}/truvari/tp-comp.vcf
#     bcftools view ${ROOT_DIR}/merged.vcf.gz > ${ROOT_DIR}/merged.vcf
#     java SupportedByCallers ${ROOT_DIR}/merged.vcf ${SVLEN_MIN} ${SVLEN_MAX} ${SVLEN_BINS} > ${ROOT_DIR}/${TR_STATUS}_callers_merged.log
#     java SupportedByCallers ${ROOT_DIR}/truvari/tp-comp.vcf ${SVLEN_MIN} ${SVLEN_MAX} ${SVLEN_BINS} > ${ROOT_DIR}/${TR_STATUS}_callers_tp.log
#     rm -f ${ROOT_DIR}/merged.vcf ${ROOT_DIR}/truvari/tp-comp.vcf
done


# # SVJEDIGRAPH
# # We regenotype the entire VCF before subsetting, to allow reads to be mapped to
# # the entire variation graph of the genome, rather than to just a subgraph.
# #
# formatVcf ${CALLER_VCF} ${ROOT_DIR}/merged.vcf.gz 0 0
# gunzip -c ${ROOT_DIR}/merged.vcf.gz > ${ROOT_DIR}/merged.vcf
# ${CLEAN_VCF_COMMAND} ${ROOT_DIR}/merged.vcf ${ROOT_DIR} ${SVLEN_MAX} 0 ${ROOT_DIR}/tmp1.vcf
# rm -f ${ROOT_DIR}/merged.vcf
# svjedi-graph.py --threads ${N_THREADS} --ref ${REFERENCE_FA} --reads ${READS_FASTQ_GZ} --vcf ${ROOT_DIR}/tmp1.vcf --prefix ${ROOT_DIR}/test
# rm -f ${ROOT_DIR}/tmp1.vcf
# bcftools view --header-only ${ROOT_DIR}/merged.vcf.gz | grep -vwE "(GT|DP|AD|PL)" > ${ROOT_DIR}/header.txt
# N_LINES=$(wc -l < ${ROOT_DIR}/header.txt)
# head -n $(( ${N_LINES} - 1 )) ${ROOT_DIR}/header.txt > ${ROOT_DIR}/tmp1.vcf
# echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> ${ROOT_DIR}/tmp1.vcf
# echo "##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Total number of informative read alignments across all alleles (after normalization for unbalanced SVs)\">" >> ${ROOT_DIR}/tmp1.vcf
# echo "##FORMAT=<ID=AD,Number=2,Type=String,Description=\"Number of informative read alignments supporting each allele (after normalization by breakpoint number for unbalanced SVs)\">" >> ${ROOT_DIR}/tmp1.vcf
# echo "##FORMAT=<ID=PL,Number=G,Type=Float,Description=\"Phred-scaled likelihood for each genotype\">" >> ${ROOT_DIR}/tmp1.vcf
# echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> ${ROOT_DIR}/tmp1.vcf
# grep '^[^#]' ${ROOT_DIR}/test_genotype.vcf >> ${ROOT_DIR}/tmp1.vcf
# bcftools sort --output-type z ${ROOT_DIR}/tmp1.vcf > ${ROOT_DIR}/regenotyped.vcf.gz
# tabix -f ${ROOT_DIR}/regenotyped.vcf.gz
# mv ${ROOT_DIR}/regenotyped.vcf.gz ${ROOT_DIR}/regenotyped-backup.vcf.gz
# mv ${ROOT_DIR}/regenotyped.vcf.gz.tbi ${ROOT_DIR}/regenotyped-backup.vcf.gz.tbi
# for TR_STATUS in 0 1 2 3; do
#     LOG_FILE="${ROOT_DIR}/${TR_STATUS}.log"
#     # Formatting truth VCF
#     truvari anno svinfo --minsize 0 ${TRUTH_VCF} | bgzip > ${ROOT_DIR}/tmp-1.vcf.gz
#     tabix ${ROOT_DIR}/tmp-1.vcf.gz
#     formatVcf ${ROOT_DIR}/tmp-1.vcf.gz ${ROOT_DIR}/truth.vcf.gz 1 ${TR_STATUS}
#     # Formatting regenotyped VCF
#     formatVcf ${ROOT_DIR}/regenotyped-backup.vcf.gz ${ROOT_DIR}/regenotyped.vcf.gz 0 ${TR_STATUS}
#     # Comparing
#     afterRegenotyping ${TR_STATUS} svjedigraph ${LOG_FILE}
# done

# # Sam's VCF
# SAMS_VCF="/Users/fcunial/Downloads/kanpig.score.xgb.vcf.gz"
# N_CALLS=$(bcftools view --no-header ${SAMS_VCF} | wc -l)
# FILTER_STRING="GT=\"./.\" || GT=\"./0\" || GT=\"0/.\" || GT=\"0/0\" || GT=\".|.\" || GT=\".|0\" || GT=\"0|.\" || GT=\"0|0\""
# N_CALLS_00_BEFORE=$(bcftools filter -i "${FILTER_STRING}" ${SAMS_VCF} | grep '^[^#]' | wc -l || echo 0)
# for TR_STATUS in 0 1 2 3; do
#     LOG_FILE="${ROOT_DIR}/${TR_STATUS}.log"
#     # Formatting truth VCF
#     truvari anno svinfo --minsize 0 ${TRUTH_VCF} | bgzip > ${ROOT_DIR}/tmp-2.vcf.gz
#     tabix -f ${ROOT_DIR}/tmp-2.vcf.gz
#     bcftools view --output-type z ${ROOT_DIR}/tmp-2.vcf.gz chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX,chrY > ${ROOT_DIR}/tmp-1.vcf.gz
#     tabix -f ${ROOT_DIR}/tmp-1.vcf.gz
#     formatVcf ${ROOT_DIR}/tmp-1.vcf.gz ${ROOT_DIR}/truth.vcf.gz 1 ${TR_STATUS}
#     # Formatting regenotyped VCF
#     formatVcf ${SAMS_VCF} ${ROOT_DIR}/regenotyped.vcf.gz 0 ${TR_STATUS}
#     # Comparing
#     afterRegenotyping ${TR_STATUS} sams ${LOG_FILE}
# done

# Outputting
# rm -f ${ROOT_DIR}/${SAMPLE_ID}*.tar.gz
# tar -czf ${ROOT_DIR}/${SAMPLE_ID}_logs.tar.gz ${ROOT_DIR}/*.log
# tar -czf ${ROOT_DIR}/${SAMPLE_ID}_zeroReads.tar.gz ${ROOT_DIR}/*_zeroReads.vcf
