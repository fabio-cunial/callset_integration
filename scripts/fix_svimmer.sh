#!/bin/bash
RAW_SVIMMER_FILE=$1

set -euxo pipefail

bcftools view --header-only ${RAW_SVIMMER_FILE} > header1.txt
N_ROWS=$(wc -l < header1.txt)
head -n $(( ${N_ROWS} - 1 )) header1.txt > fixed.vcf
echo "##FILTER=<ID=GT,Description="Genotype filter">" >> fixed.vcf
echo "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">" >> fixed.vcf
echo "##FORMAT=<ID=DR,Number=1,Type=Integer,Description="Number of reference reads">" >> fixed.vcf
echo "##FORMAT=<ID=DV,Number=1,Type=Integer,Description="Number of variant reads">" >> fixed.vcf
echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> fixed.vcf

bcftools view --no-header ${RAW_SVIMMER_FILE} > body1.txt
cat body1.txt | awk '{
	format="GT:GQ:DR:DV"; \
	n_samples=1; \
	printf("%s",$1); \
	for (i=2; i<=NF; i++) printf("\t%s",$i); \
	printf("\t%s",format); \
	for (i=1; i<=n_samples; i++) printf("\t0/1:99:1:99"); \
	printf("\n"); \
}' >> fixed.vcf