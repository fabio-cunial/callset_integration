version 1.0


#
workflow GetRegenotypedTruePositives {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File truth_bed
        File reference_fa
        File reference_fai
        Int svlen_min
        Int svlen_max
    }
    parameter_meta {
    }
    
    call GetRegenotypedTruePositivesImpl {
        input:
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            truth_vcf_gz = truth_vcf_gz,
            truth_vcf_gz_tbi = truth_vcf_gz_tbi,
            truth_bed = truth_bed,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            svlen_min = svlen_min,
            svlen_max = svlen_max
    }
    
    output {
        File regenotyped_sniffles = GetRegenotypedTruePositivesImpl.regenotyped_sniffles
        File regenotyped_sniffles_tbi = GetRegenotypedTruePositivesImpl.regenotyped_sniffles_tbi
        File regenotyped_sniffles_tp_comp = GetRegenotypedTruePositivesImpl.regenotyped_sniffles_tp_comp
        File regenotyped_sniffles_tp_comp_tbi = GetRegenotypedTruePositivesImpl.regenotyped_sniffles_tp_comp_tbi
        File regenotyped_sniffles_tp_base = GetRegenotypedTruePositivesImpl.regenotyped_sniffles_tp_base
        File regenotyped_sniffles_tp_base_tbi = GetRegenotypedTruePositivesImpl.regenotyped_sniffles_tp_base_tbi
        File regenotyped_kanpig = GetRegenotypedTruePositivesImpl.regenotyped_kanpig
        File regenotyped_kanpig_tbi = GetRegenotypedTruePositivesImpl.regenotyped_kanpig_tbi
        File regenotyped_kanpig_tp_comp = GetRegenotypedTruePositivesImpl.regenotyped_kanpig_tp_comp
        File regenotyped_kanpig_tp_comp_tbi = GetRegenotypedTruePositivesImpl.regenotyped_kanpig_tp_comp_tbi
        File regenotyped_kanpig_tp_base = GetRegenotypedTruePositivesImpl.regenotyped_kanpig_tp_base
        File regenotyped_kanpig_tp_base_tbi = GetRegenotypedTruePositivesImpl.regenotyped_kanpig_tp_base_tbi
    }
}


task GetRegenotypedTruePositivesImpl {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File truth_bed
        File reference_fa
        File reference_fai
        Int svlen_min
        Int svlen_max
    }
    parameter_meta {
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
        KANPIG_SIZEMAX="10000"  # From Adam's suggestion:
        # "I noticed that you used for kanpig --sizemax 1000000 . You're going to get lower recall with that. Currently kanpig is using a very naive clustering strategy to figure out which variants should be considered together. The boundaries of the variant graphs are set to min_start/max_end and pileups are made over the region. Since kanpig is also only looking at pileups of reads that span the region, large variants can preclude smaller variants from getting a chance to have read support. I'm working on better clustering and not needing only spanning reads, but for now sizemax should be set to something like 10k, or maybe even ~75% of the mean insert size."
        chmod +x ~{docker_dir}/kanpig

        # Makes sure that the truth and merged VCF are in the right format and
        # contain only a specific set of calls.
        #
        function formatVcf() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF_GZ=$2
            local IS_TRUTH=$3
            
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

            # If not a truth VCF: keeping only calls in the given length range.
            if [ ${IS_TRUTH} -eq 1 ]; then
                cp tmp2.vcf.gz tmp3.vcf.gz
                cp tmp2.vcf.gz.tbi tmp3.vcf.gz.tbi
            else
                rm -rf tmp3.vcf*
                FILTER_STRING="SVLEN>=~{svlen_min} && SVLEN<=~{svlen_max}"
                bcftools filter -i "${FILTER_STRING}" --output-type z tmp2.vcf.gz > tmp3.vcf.gz
                tabix -f tmp3.vcf.gz
            fi

            # Finalizing
            cp tmp3.vcf.gz ${OUTPUT_VCF_GZ}
            cp tmp3.vcf.gz.tbi ${OUTPUT_VCF_GZ}.tbi
            rm -f tmp*.vcf*
        }


        # Main program
        rm -f ~{alignments_bai}
        samtools index -@ ${N_THREADS} ~{alignments_bam}
    
        # Formatting truth and merged VCF
        formatVcf ~{truvari_collapsed_vcf_gz} merged.vcf.gz 0
        truvari anno svinfo --minsize 0 ~{truth_vcf_gz} | bgzip > tmp-1.vcf.gz
        tabix -f tmp-1.vcf.gz
        formatVcf tmp-1.vcf.gz truth.vcf.gz 1

        # SNIFFLES FORCE
        rm -f regenotyped.vcf.gz
        sniffles --threads ${N_THREADS} --reference ~{reference_fa} --input ~{alignments_bam} --genotype-vcf merged.vcf.gz --vcf regenotyped_sniffles.vcf.gz
        tabix -f regenotyped_sniffles.vcf.gz
        rm -rf truvari/
        truvari bench ${TRUVARI_BENCH_SETTINGS} -b truth.vcf.gz -c regenotyped_sniffles.vcf.gz -o truvari/
        bcftools sort --output-type z truvari/tp-comp.vcf.gz > regenotyped_sniffles_tp_comp.vcf.gz
        tabix -f regenotyped_sniffles_tp_comp.vcf.gz
        bcftools sort --output-type z truvari/tp-base.vcf.gz > regenotyped_sniffles_tp_base.vcf.gz
        tabix -f regenotyped_sniffles_tp_base.vcf.gz

        # KANPIG
        export RUST_BACKTRACE=1
        ~{docker_dir}/kanpig --threads ${N_THREADS} --sizemin ~{svlen_min} --sizemax ${KANPIG_SIZEMAX} --input merged.vcf.gz --bam ~{alignments_bam} --reference ~{reference_fa} --out tmp1.vcf.gz
        bcftools sort --output-type z tmp1.vcf.gz > regenotyped_kanpig.vcf.gz
        tabix -f regenotyped_kanpig.vcf.gz
        rm -f tmp1.vcf.gz
        rm -rf truvari/
        truvari bench ${TRUVARI_BENCH_SETTINGS} -b truth.vcf.gz -c regenotyped_kanpig.vcf.gz -o truvari/
        bcftools sort --output-type z truvari/tp-comp.vcf.gz > regenotyped_kanpig_tp_comp.vcf.gz
        tabix -f regenotyped_kanpig_tp_comp.vcf.gz
        bcftools sort --output-type z truvari/tp-base.vcf.gz > regenotyped_kanpig_tp_base.vcf.gz
        tabix -f regenotyped_kanpig_tp_base.vcf.gz
    >>>

    output {
        File regenotyped_sniffles = work_dir + "/regenotyped_sniffles.vcf.gz"
        File regenotyped_sniffles_tbi = work_dir + "/regenotyped_sniffles.vcf.gz.tbi"
        File regenotyped_sniffles_tp_comp = work_dir + "/regenotyped_sniffles_tp_comp.vcf.gz"
        File regenotyped_sniffles_tp_comp_tbi = work_dir + "/regenotyped_sniffles_tp_comp.vcf.gz.tbi"
        File regenotyped_sniffles_tp_base = work_dir + "/regenotyped_sniffles_tp_base.vcf.gz"
        File regenotyped_sniffles_tp_base_tbi = work_dir + "/regenotyped_sniffles_tp_base.vcf.gz.tbi"
        File regenotyped_kanpig = work_dir + "/regenotyped_kanpig.vcf.gz"
        File regenotyped_kanpig_tbi = work_dir + "/regenotyped_kanpig.vcf.gz.tbi"
        File regenotyped_kanpig_tp_comp = work_dir + "/regenotyped_kanpig_tp_comp.vcf.gz"
        File regenotyped_kanpig_tp_comp_tbi = work_dir + "/regenotyped_kanpig_tp_comp.vcf.gz.tbi"
        File regenotyped_kanpig_tp_base = work_dir + "/regenotyped_kanpig_tp_base.vcf.gz"
        File regenotyped_kanpig_tp_base_tbi = work_dir + "/regenotyped_kanpig_tp_base.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "32GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
