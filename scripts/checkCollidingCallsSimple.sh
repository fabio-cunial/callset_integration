#!/bin/bash
#
INPUT_VCF=$1
INPUT_COORDINATE=$2
N_SAMPLES="1074"


bcftools view --no-header ${INPUT_VCF} ${INPUT_COORDINATE} > tmp.vcf
for i in $(seq 1 ${N_SAMPLES} ); do
	cut -f $(( 9 + ${i} )) tmp.vcf > haps.txt
	N_01=$(grep "0|1" haps.txt | wc -l)
	if [ ${N_01} -gt 1 ]; then
		echo "${N_01} 0|1 on the ${i}-th sample:"
		cut -f 1,2,3,4,5,6,7,8,9,$(( 9 + ${i} )) tmp.vcf
        exit
	fi
	N_10=$(grep "1|0" haps.txt | wc -l)
	if [ ${N_10} -gt 1 ]; then
		echo "${N_10} 1|0 conflicts on the ${i}-th sample:"
		cut -f 1,2,3,4,5,6,7,8,9,$(( 9 + ${i} )) tmp.vcf
        exit
	fi
	N_11=$(grep "1|1" haps.txt | wc -l)
	if [ ${N_11} -gt 1 ]; then
		echo "${N_11} 1|1 conflicts on the ${i}-th sample:"
		cut -f 1,2,3,4,5,6,7,8,9,$(( 9 + ${i} )) tmp.vcf
        exit
	fi
    N_01_11=$(( ${N_01} + ${N_11} ))
    if [ ${N_01_11} -gt 1 ]; then
		echo "${N_01_11} 0|1+1|1 conflicts on the ${i}-th sample:"
		cut -f 1,2,3,4,5,6,7,8,9,$(( 9 + ${i} )) tmp.vcf
        exit
    fi
    N_10_11=$(( ${N_10} + ${N_11} ))
    if [ ${N_10_11} -gt 1 ]; then
		echo "${N_10_11} 1|0+1|1 conflicts on the ${i}-th sample:"
		cut -f 1,2,3,4,5,6,7,8,9,$(( 9 + ${i} )) tmp.vcf
        exit
    fi
    rm -f haps.txt
done
#rm -f tmp.vcf
