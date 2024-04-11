version 1.0


#
workflow GetRegenotypedVcfKanpig {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int svlen_min
        Int svlen_max
    }
    parameter_meta {
    }
    
    call GetRegenotypedVcfImpl {
        input:
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            svlen_min = svlen_min,
            svlen_max = svlen_max
    }
    
    output {
        File regenotyped_kanpig = GetRegenotypedVcfImpl.regenotyped_kanpig
        File regenotyped_kanpig_tbi = GetRegenotypedVcfImpl.regenotyped_kanpig_tbi
    }
}


task GetRegenotypedVcfImpl {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int svlen_min
        Int svlen_max
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int mem_gb = 60
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{mem_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        # "I noticed that you used for kanpig --sizemax 1000000 . You're going to get lower recall with that. Currently kanpig is using a very naive clustering strategy to figure out which variants should be considered together. The boundaries of the variant graphs are set to min_start/max_end and pileups are made over the region. Since kanpig is also only looking at pileups of reads that span the region, large variants can preclude smaller variants from getting a chance to have read support. I'm working on better clustering and not needing only spanning reads, but for now sizemax should be set to something like 10k, or maybe even ~75% of the mean insert size."
        KANPIG_SIZEMAX="10000"  # From Adam's suggestion above
        # I tried a dozen different conditions and these params still seem to
        # have the best balance on 8x and a collapsed single-sample vcf:
        KANPIG_PARAMS="--chunksize 1000 --sizesim 0.90 --seqsim 0.85 --hapsim 0.9999 --no-prune --maxpaths 10000"
        chmod +x ~{docker_dir}/kanpig

        # Makes sure that the merged VCF is in the right format and contains
        # only a specific set of calls.
        #
        function formatVcf() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF_GZ=$2
            
            # Removing multiallelic records
            rm -f tmp0.vcf.gz*
            bcftools norm --multiallelics - --output-type z ${INPUT_VCF_GZ} > tmp0.vcf.gz
            tabix -f tmp0.vcf.gz
            
            # Adding REF and ALT with Adam's script, and removing BNDs.
            rm -f tmp1.vcf*
            python ~{docker_dir}/resolve.py tmp0.vcf.gz ~{reference_fa} | bcftools view -i "SVTYPE!='BND'" | bcftools sort -O z -o tmp1.vcf.gz
            tabix -f tmp1.vcf.gz

            # Keeping only calls in the given length range.
            rm -rf tmp2.vcf*
            FILTER_STRING="SVLEN>=~{svlen_min} && SVLEN<=~{svlen_max}"
            bcftools filter -i "${FILTER_STRING}" --output-type z tmp1.vcf.gz > tmp2.vcf.gz
            tabix -f tmp2.vcf.gz

            # Finalizing
            cp tmp2.vcf.gz ${OUTPUT_VCF_GZ}
            cp tmp2.vcf.gz.tbi ${OUTPUT_VCF_GZ}.tbi
            rm -f tmp*.vcf*
        }


        # Main program
        #rm -f ~{alignments_bai}
        #samtools index -@ ${N_THREADS} ~{alignments_bam}
    
        # Formatting the merged VCF
        formatVcf ~{truvari_collapsed_vcf_gz} merged.vcf.gz 0

        # KANPIG
        export RUST_BACKTRACE=1
        ~{docker_dir}/kanpig --threads ${N_THREADS} --sizemin ~{svlen_min} --sizemax ~{svlen_max} ${KANPIG_PARAMS} --input merged.vcf.gz --bam ~{alignments_bam} --reference ~{reference_fa} --out tmp1.vcf.gz
        bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z tmp1.vcf.gz > regenotyped_kanpig.vcf.gz
        tabix -f regenotyped_kanpig.vcf.gz
        rm -f tmp1.vcf.gz
    >>>

    output {
        File regenotyped_kanpig = work_dir + "/regenotyped_kanpig.vcf.gz"
        File regenotyped_kanpig_tbi = work_dir + "/regenotyped_kanpig.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: mem_gb + "GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
