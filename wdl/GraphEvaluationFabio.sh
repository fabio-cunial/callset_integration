#!/bin/bash
#
# -------------------------------- INPUT ---------------------------------------
CHROMOSOME=$1
INTERVAL_MAX_LENGTH="30000"
FLANK_LENGTH="150"
EVALUATION_BEDS=$(ls /mnt/disks/disk-1/experiments/evaluation_beds/*.bed)
VCF_LIST_FILE="/mnt/disks/disk-1/experiments/vcf_list.txt"
REFERENCE_FA="/mnt/disks/disk-1/experiments/reference/chm13v2.0.ebv.fa"
TANDEMS_BED="/mnt/disks/disk-1/experiments/reference/human_chm13v2.0_maskedY_rCRS.trf.bed"
BAM_CSV="/mnt/disks/disk-1/experiments/hprc_94_haps_vs_chm13_v2.csv"
HAPESTRY_DIR="/home/fcunial/sv_merge/build"
TIME_COMMAND="time --verbose"
# ------------------------------------------------------------------------------


set -euxo pipefail

N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))


# Downloading all the calls in $chromosome$.
export GCS_OAUTH_TOKEN=$(gcloud auth print-access-token)
VCFS=""; CLUSTER_BY=""; TOOLS=""; LINE="0";
while read ROW; do
    TOOL_ID=${ROW%,*}
    TOOL_FILE=${ROW#*,}
    bcftools view --output-type z ${TOOL_FILE} ${CHROMOSOME} > tmp.vcf.gz
    tabix -f tmp.vcf.gz
    bcftools norm --multiallelics - --output-type v tmp.vcf.gz > ${TOOL_ID}.vcf
    rm -f tmp.vcf.gz
    LINE=$(( ${LINE} + 1 ))
    if [[ ${LINE} -eq 1 ]]; then
        VCFS=${TOOL_ID}.vcf
        CLUSTER_BY=${TOOL_ID}.vcf
        TOOLS=${TOOL_ID}
    else
        VCFS="${VCFS},${TOOL_ID}.vcf"
        TOOLS="${TOOLS} ${TOOL_ID}"
    fi
done < ${VCF_LIST_FILE}

# Evaluating
EVALUATION_NAME="${CHROMOSOME}_evaluation"
ANALYSIS_NAME="${CHROMOSOME}_analysis"
rm -rf ./${EVALUATION_NAME}
${TIME_COMMAND} ${HAPESTRY_DIR}/evaluate \
    --n_threads ${N_THREADS} \
    --output_dir ./${EVALUATION_NAME} \
    --bam_csv ${BAM_CSV} \
    --vcfs ${VCFS} \
    --cluster_by ${CLUSTER_BY} \
    --tandems ${TANDEMS_BED} \
    --ref ${REFERENCE_FA} \
    --interval_max_length ${INTERVAL_MAX_LENGTH} \
    --flank_length ${FLANK_LENGTH} \
    --debug
TOOLS=$(echo "${TOOLS}" | tr '.' '_')
rm -rf ./${ANALYSIS_NAME}
${TIME_COMMAND} ${HAPESTRY_DIR}/analyze_evaluation \
    --input_dir ./${EVALUATION_NAME} \
    --output_dir ./${ANALYSIS_NAME} \
    --tools ${TOOLS} \
    --beds ${EVALUATION_BEDS}

export GZIP=-1
${TIME_COMMAND} tar -czf ${ANALYSIS_NAME}.tar.gz ./${ANALYSIS_NAME}
${TIME_COMMAND} tar -czf ${EVALUATION_NAME}.tar.gz --exclude='*.fasta' --exclude='*.fa' ./${EVALUATION_NAME}
