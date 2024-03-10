#!/bin/bash
#
# Removes duplicated remote BAM filenames.
# Assumes that every ID in the input CSV is distinct.
#
set -euxo pipefail

INPUT_CSV="hgsvc2_28_haps_vc_grch38.csv"
OUTPUT_CSV="hgsvc2_28_haps_vc_grch38_new.csv"

rm -f ${OUTPUT_CSV}
while read ROW; do
    ID=${ROW%,*}
    ADDRESS=${ROW#*,}
    PREFIX=${ADDRESS%/*}
    gsutil mv ${ADDRESS} ${PREFIX}/${ID}.bam
    gsutil mv ${ADDRESS}.bai ${PREFIX}/${ID}.bam.bai
    echo ${ID},${PREFIX}/${ID}.bam >> ${OUTPUT_CSV}
done < ${INPUT_CSV}