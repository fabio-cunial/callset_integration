#!/bin/bash
#
WINDOW=$1
RUN_ID=$2

RUNS_ROOT_DIR="/mnt/disks/disk-1/local_hapestry_run"
VCF_GZ="/mnt/disks/disk-1/local_hapestry_run/hprc_chm13.vcf.gz"
BAM_CSV="/mnt/disks/disk-1/local_hapestry_run/hprc_47_8x_reads.csv"
TANDEMS_BED="/mnt/disks/disk-1/local_hapestry_run/human_chm13v2.0_maskedY_rCRS.trf.bed"
REFERENCE_FA="/mnt/disks/disk-1/local_hapestry_run/chm13v2.0.ebv.fa"
INTERVAL_MAX_LENGTH="15000"
FLANK_LENGTH="200"
MIN_SV_LENGTH="1"
GRAPHALIGNER_TIMEOUT="400"
SOLVER_TIMEOUT="900"
MIN_READ_HAP_IDENTITY="0"
D_WEIGHT="32"
N_THREADS="1"
HAPESTRY_COMMAND="/home/fcunial/sv_merge/build/hapestry"
TIME_COMMAND="/usr/bin/time --verbose"

RUN_DIR=${RUNS_ROOT_DIR}/${RUN_ID}
rm -f tmp.vcf
bcftools view --output-type v ${VCF_GZ} ${WINDOW} > tmp.vcf
rm -rf ${RUN_DIR}
${TIME_COMMAND} ${HAPESTRY_COMMAND} --output_dir ${RUN_DIR} \
        --bam_csv ${BAM_CSV} \
        --vcf tmp.vcf \
        --tandems ${TANDEMS_BED} \
        --ref ${REFERENCE_FA} \
        --interval_max_length ${INTERVAL_MAX_LENGTH} \
        --min_sv_length ${MIN_SV_LENGTH} \
        --flank_length ${FLANK_LENGTH} \
        --graphaligner_timeout ${GRAPHALIGNER_TIMEOUT} \
        --solver_timeout ${SOLVER_TIMEOUT} \
        --min_read_hap_identity ${MIN_READ_HAP_IDENTITY} \
        --d_weight ${D_WEIGHT} \
        --n_threads ${N_THREADS} \
        --bam_not_hardclipped
