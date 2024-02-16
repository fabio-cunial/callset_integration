#!/bin/bash
#

for CHROMOSOME in chrY chrX chr22 chr21 chr20 chr19 chr18 ; do
    ./GraphEvaluationFabio.sh ${CHROMOSOME}
done
