#!/bin/bash
#
ROOT_DIR="gs://fc-bb12940d-f7ba-4515-8a98-42de84f63c34/submissions/590a3fc4-fdb5-485a-8bfd-f345b128a6ad/GraphEvaluation/5bff8350-addb-43f2-acb8-c1201aa64761/call-EvaluateChromosome"
for i in $(seq 21 -1 5); do
    gsutil -m cp ${ROOT_DIR}/shard-$(( i - 1 ))/hapestry/chr${i}_evaluation.tar.gz .
    tar -xzf chr${i}_evaluation.tar.gz
    rm -f chr${i}_evaluation.tar.gz
    cp ./chr${i}_evaluation/windows_omitted.bed windows_omitted_${i}.bed 
    rm -rf ./chr${i}_evaluation/
done
