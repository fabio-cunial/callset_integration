version 1.0


#
workflow ROC {
    input {
        String sample_id
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reads_fastq_gz
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File truth_bed
        File reference_fa
        File reference_fai
        File tr_bed
        File nontr_bed
        Float tr_small_overlap
        Float tr_large_overlap
        String svlen_bins 
        Int svlen_min
        Int svlen_max
        File chromosome_files_gz
    }
    parameter_meta {
        svlen_bins: "Comma-separated, sorted."
    }
    
    call ROCImpl {
        input:
            sample_id = sample_id,
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reads_fastq_gz = reads_fastq_gz,
            truth_vcf_gz = truth_vcf_gz,
            truth_vcf_gz_tbi = truth_vcf_gz_tbi,
            truth_bed = truth_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            tr_bed = tr_bed,
            nontr_bed = nontr_bed,
            tr_small_overlap = tr_small_overlap,
            tr_large_overlap = tr_large_overlap,
            svlen_bins = svlen_bins,
            svlen_min = svlen_min,
            svlen_max = svlen_max,
            chromosome_files_gz = chromosome_files_gz
    }
    
    output {
        File logs = ROCImpl.logs
        File zero_reads = ROCImpl.zero_reads
    }
}


task ROCImpl {
    input {
        String sample_id
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reads_fastq_gz
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File truth_bed
        File reference_fa
        File reference_fai
        File tr_bed
        File nontr_bed
        Float tr_small_overlap
        Float tr_large_overlap
        String svlen_bins 
        Int svlen_min
        Int svlen_max
        File chromosome_files_gz
    }
    parameter_meta {
        svlen_bins: "Comma-separated, sorted."
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TRUVARI_BENCH_SETTINGS="--sizefilt ~{svlen_min} --sizemax ~{svlen_max} --sizemin 0"
        chmod +x ~{docker_dir}/kanpig

        # Makes sure that the truth and merged VCF are in the right format and
        # contain only a specific set of calls.
        #
        # @param 4: keep only calls with: 0=no filter; 1=large overlap with the
        # TR track; 2=small overlap with the TR track; 3=large overlap with the
        # non-TR track.
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
            rm -f tmp1.vcf*
            python ~{docker_dir}/resolve.py tmp0.vcf.gz ~{reference_fa} | bcftools view -i "SVTYPE!='BND'" | bcftools sort -O z -o tmp1.vcf.gz
            tabix -f tmp1.vcf.gz

            # Keeping only calls in the truth BED (any overlap).
            rm -rf tmp2.vcf*
            bcftools view --header-only tmp1.vcf.gz > tmp2.vcf
            bedtools intersect -u -a tmp1.vcf.gz -b ~{truth_bed} >> tmp2.vcf
            rm -f tmp2.vcf.gz*; bgzip tmp2.vcf
            tabix -f tmp2.vcf.gz

            # If not a truth VCF: keeping only calls in the given TR track
            # (with min. overlap).
            if [ ${IS_TRUTH} -eq 1 ]; then
                cp tmp2.vcf.gz tmp3.vcf.gz
                cp tmp2.vcf.gz.tbi tmp3.vcf.gz.tbi
            else
                rm -rf tmp3.vcf*
                if [ ${TR_STATUS} -eq 0 ]; then
                    cp tmp2.vcf.gz tmp3.vcf.gz
                    cp tmp2.vcf.gz.tbi tmp3.vcf.gz.tbi
                elif [ ${TR_STATUS} -eq 1 ]; then
                    bcftools view --header-only tmp2.vcf.gz > tmp3.vcf
                    bedtools intersect -u -f ~{tr_large_overlap} -a tmp2.vcf.gz -b ~{tr_bed} >> tmp3.vcf
                    rm -f tmp3.vcf.gz*; bgzip tmp3.vcf
                    tabix -f tmp3.vcf.gz
                elif [ ${TR_STATUS} -eq 2 ]; then
                    bcftools view --header-only tmp2.vcf.gz > tmp3.vcf
                    bedtools intersect -u -f ~{tr_small_overlap} -a tmp2.vcf.gz -b ~{tr_bed} >> tmp3.vcf
                    rm -f tmp3.vcf.gz*; bgzip tmp3.vcf
                    tabix -f tmp3.vcf.gz
                elif [ ${TR_STATUS} -eq 3 ]; then
                    bcftools view --header-only tmp2.vcf.gz > tmp3.vcf
                    bedtools intersect -u -f ~{tr_large_overlap} -a tmp2.vcf.gz -b ~{nontr_bed} >> tmp3.vcf
                    rm -f tmp3.vcf.gz*; bgzip tmp3.vcf
                    tabix -f tmp3.vcf.gz
                else
                    :
                fi
            fi

            # If not a truth VCF: keeping only calls in the given length range.
            if [ ${IS_TRUTH} -eq 1 ]; then
                cp tmp3.vcf.gz tmp4.vcf.gz
                cp tmp3.vcf.gz.tbi tmp4.vcf.gz.tbi
            else
                rm -rf tmp4.vcf*
                FILTER_STRING="SVLEN>=~{svlen_min} && SVLEN<=~{svlen_max}"
                bcftools filter -i "${FILTER_STRING}" --output-type z tmp3.vcf.gz > tmp4.vcf.gz
                tabix tmp4.vcf.gz
            fi

            # Finalizing
            cp tmp4.vcf.gz ${OUTPUT_VCF_GZ}
            cp tmp4.vcf.gz.tbi ${OUTPUT_VCF_GZ}.tbi
            rm -f tmp*.vcf*
        }


        # @param 3: 0=no filter; 1=large overlap with the TR track; 2=small
        # overlap with the TR track; 3=large overlap with the non-TR track.
        #
        function afterRegenotyping() {
            local TR_STATUS=$1
            local GENOTYPER=$2
            local LOG_FILE=$3
    
            tabix -f regenotyped.vcf.gz
            N_CALLS_00_AFTER=$(bcftools filter -i "${FILTER_STRING}" regenotyped.vcf.gz | grep '^[^#]' | wc -l)
            VALUE_BEFORE="${N_CALLS_00_BEFORE}/${N_CALLS}"; VALUE_BEFORE=$(bc -l <<< ${VALUE_BEFORE})
            VALUE_AFTER="${N_CALLS_00_AFTER}/${N_CALLS}"; VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
            echo "${TR_STATUS},${GENOTYPER} Fraction of 0/0 calls: Before regenotyping: ${VALUE_BEFORE} -> After regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
            if [ ${GENOTYPER} = sniffles -o ${GENOTYPER} = cutesv ]; then
                N_CALLS_DV_ZERO_AFTER=$(bcftools filter -i 'DV=0' regenotyped.vcf.gz | grep '^[^#]' | wc -l)
                VALUE_AFTER="${N_CALLS_DV_ZERO_AFTER}/${N_CALLS}"; VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
                echo "${TR_STATUS},${GENOTYPER} Fraction of calls with 0 supporting reads after regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
                bcftools view merged.vcf.gz > merged.vcf
                bcftools view regenotyped.vcf.gz > regenotyped.vcf
                java -cp ~{docker_dir} SupportedByZeroReads DV merged.vcf regenotyped.vcf ~{svlen_min} ~{svlen_max} ~{svlen_bins} > ${TR_STATUS}_${GENOTYPER}_zeroReads.log
                rm -f merged.vcf regenotyped.vcf
            elif [ ${GENOTYPER} = kanpig ]; then
                N_CALLS_DV_ZERO_AFTER=$(bcftools filter -i 'FORMAT/AD[0:1]=0' regenotyped.vcf.gz | grep '^[^#]' | wc -l)
                VALUE_AFTER="${N_CALLS_DV_ZERO_AFTER}/${N_CALLS}"; VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
                echo "${TR_STATUS},${GENOTYPER} Fraction of calls with 0 supporting reads after regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
                bcftools view merged.vcf.gz > merged.vcf
                bcftools view regenotyped.vcf.gz > regenotyped.vcf
                java -cp ~{docker_dir} SupportedByZeroReads AD merged.vcf regenotyped.vcf ~{svlen_min} ~{svlen_max} ~{svlen_bins} > ${TR_STATUS}_${GENOTYPER}_zeroReads.log
                rm -f merged.vcf regenotyped.vcf
            elif [ ${GENOTYPER} = svjedigraph ]; then
                N_CALLS_DV_ZERO_AFTER=$(bcftools filter -i 'FORMAT/AD[0:1]="0"' regenotyped.vcf.gz | grep '^[^#]' | wc -l)
                VALUE_AFTER="${N_CALLS_DV_ZERO_AFTER}/${N_CALLS}"; VALUE_AFTER=$(bc -l <<< ${VALUE_AFTER})
                echo "${TR_STATUS},${GENOTYPER} Fraction of calls with 0 supporting reads after regenotyping: ${VALUE_AFTER}" >> ${LOG_FILE}
                bcftools view merged.vcf.gz > merged.vcf
                bcftools view regenotyped.vcf.gz > regenotyped.vcf
                java -cp ~{docker_dir} SupportedByZeroReads AD merged.vcf regenotyped.vcf ~{svlen_min} ~{svlen_max} ~{svlen_bins} > ${TR_STATUS}_${GENOTYPER}_zeroReads.log
                rm -f merged.vcf regenotyped.vcf
            fi

            # Computing TPs and ROC curves
            rm -rf truvari/
            truvari bench ${TRUVARI_BENCH_SETTINGS} -b truth.vcf.gz -c regenotyped.vcf.gz -o truvari/
            bcftools view regenotyped.vcf.gz > regenotyped.vcf
            bcftools view truvari/tp-comp.vcf.gz > truvari/tp-comp.vcf
            bcftools view truth.vcf.gz > truth.vcf
            if [ ${GENOTYPER} = sniffles -o ${GENOTYPER} = cutesv ]; then
                java -cp ~{docker_dir} ROCcurve GQ 1 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve GQ 0 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve DV 1 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve DV 0 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                bcftools filter -i 'DV=0' regenotyped.vcf.gz > ${TR_STATUS}_${GENOTYPER}_zeroReads.vcf
            elif [ ${GENOTYPER} = kanpig ]; then
                java -cp ~{docker_dir} ROCcurve SQ 1 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve SQ 0 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve AD 1 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve AD 0 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                bcftools filter -i 'FORMAT/AD[0:1]=0' regenotyped.vcf.gz > ${TR_STATUS}_${GENOTYPER}_zeroReads.vcf
            elif [ ${GENOTYPER} = svjedigraph ]; then
                #java -cp ~{docker_dir} ROCcurve PL 1 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                #java -cp ~{docker_dir} ROCcurve PL 0 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve AD 1 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                java -cp ~{docker_dir} ROCcurve AD 0 ~{svlen_min} ~{svlen_max} regenotyped.vcf truvari/tp-comp.vcf truth.vcf ~{work_dir} ${TR_STATUS} ${GENOTYPER} ~{svlen_bins}
                bcftools filter -i 'FORMAT/AD[0:1]="0"' regenotyped.vcf.gz > ${TR_STATUS}_${GENOTYPER}_zeroReads.vcf
            fi
        }


        # Main program
        tar -xzf ~{chromosome_files_gz} -C .
        for TR_STATUS in 0 1 2 3; do
            LOG_FILE="${TR_STATUS}.log"
            rm -f ${LOG_FILE}
    
            # Formatting truth and merged VCF
            formatVcf ~{truvari_collapsed_vcf_gz} merged.vcf.gz 0 ${TR_STATUS}
            truvari anno svinfo --minsize 0 ~{truth_vcf_gz} | bgzip > tmp0.vcf.gz
            tabix tmp0.vcf.gz
            formatVcf tmp0.vcf.gz truth.vcf.gz 1 ${TR_STATUS}

            # Regenotyping
            rm -f regenotyped.vcf.gz
            N_CALLS=$(bcftools view --no-header merged.vcf.gz | wc -l)
            FILTER_STRING="GT=\"./.\" || GT=\"./0\" || GT=\"0/.\" || GT=\"0/0\" || GT=\".|.\" || GT=\".|0\" || GT=\"0|.\" || GT=\"0|0\""
            N_CALLS_00_BEFORE=$(bcftools filter -i "${FILTER_STRING}" merged.vcf.gz | grep '^[^#]' | wc -l)

            # SNIFFLES FORCE
            sniffles --threads ${N_THREADS} --reference ~{reference_fa} --input ~{alignments_bam} --genotype-vcf merged.vcf.gz --vcf regenotyped.vcf.gz
            afterRegenotyping ${TR_STATUS} sniffles ${LOG_FILE}

            # CUTESV FORCE
            rm -rf ./cutesv_tmp ; mkdir ./cutesv_tmp
            cuteSV --threads ${N_THREADS} -Ivcf merged.vcf.gz --genotype --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --min_size 1 --max_size -1 ~{alignments_bam} ~{reference_fa} tmp1.vcf.gz ./cutesv_tmp
            bcftools sort --output-type z tmp1.vcf.gz > regenotyped.vcf.gz
            rm -f tmp1.vcf.gz
            afterRegenotyping ${TR_STATUS} cutesv ${LOG_FILE}

            # KANPIG
            ~{docker_dir}/kanpig --threads ${N_THREADS} --sizemin ~{svlen_min} --sizemax ~{svlen_max} --input merged.vcf.gz --bam ~{alignments_bam} --reference ~{reference_fa} --out tmp1.vcf.gz
            bcftools sort --output-type z tmp1.vcf.gz > regenotyped.vcf.gz
            rm -f tmp1.vcf.gz
            afterRegenotyping ${TR_STATUS} kanpig ${LOG_FILE}
    
            # LRCALLER
            #
            #

            # SVJEDIGRAPH
            #gunzip -c merged.vcf.gz > merged.vcf
            #java -cp ~{docker_dir} CleanVCF merged.vcf . ~{svlen_max} 0 tmp1.vcf
            #rm -f merged.vcf
            #svjedi-graph.py --threads ${N_THREADS} --ref ~{reference_fa} --reads ~{reads_fastq_gz} --vcf tmp1.vcf --prefix test
            #rm -f tmp1.vcf
            #bcftools view --header-only merged.vcf.gz | grep -vwE "(GT|DP|AD|PL)" > header.txt
            #N_LINES=$(wc -l < header.txt)
            #head -n $(( ${N_LINES} - 1 )) header.txt > tmp1.vcf
            #echo "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">" >> tmp1.vcf
            #echo "##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Total number of informative read alignments across all alleles (after normalization for unbalanced SVs)\">" >> tmp1.vcf
            #echo "##FORMAT=<ID=AD,Number=2,Type=String,Description=\"Number of informative read alignments supporting each allele (after normalization by breakpoint number for unbalanced SVs)\">" >> tmp1.vcf
            #echo "##FORMAT=<ID=PL,Number=G,Type=Float,Description=\"Phred-scaled likelihood for each genotype\">" >> tmp1.vcf
            #echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> tmp1.vcf
            #grep '^[^#]' test_genotype.vcf >> tmp1.vcf
            #bcftools sort --output-type z tmp1.vcf > regenotyped.vcf.gz
            #afterRegenotyping ${TR_STATUS} svjedigraph ${LOG_FILE}

            # Callers supporting a call
            rm -rf truvari/
            truvari bench ${TRUVARI_BENCH_SETTINGS} -b truth.vcf.gz -c merged.vcf.gz -o truvari/
            bcftools view truvari/tp-comp.vcf.gz > truvari/tp-comp.vcf
            bcftools view merged.vcf.gz > merged.vcf
            java -cp ~{docker_dir} SupportedByCallers merged.vcf ~{svlen_min} ~{svlen_max} ~{svlen_bins} > ${TR_STATUS}_callers_merged.log
            java -cp ~{docker_dir} SupportedByCallers truvari/tp-comp.vcf ~{svlen_min} ~{svlen_max} ~{svlen_bins} > ${TR_STATUS}_callers_tp.log
            rm -f merged.vcf truvari/tp-comp.vcf
        done

        # Outputting
        tar -czf ~{sample_id}_logs.tar.gz *.log
        tar -czf ~{sample_id}_zeroReads.tar.gz *_zeroReads.vcf
    >>>

    output {
        File logs = work_dir + "/" + sample_id + "_logs.tar.gz"
        File zero_reads = work_dir + "/" + sample_id + "_zeroReads.tar.gz"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "32GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
