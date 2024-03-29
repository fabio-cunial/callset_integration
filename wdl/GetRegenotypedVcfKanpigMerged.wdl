version 1.0


#
workflow GetRegenotypedVcfKanpigMerged {
    input {
        String sample_id
        File merged_vcf_gz
        File merged_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    call GetRegenotypedVcfImpl {
        input:
            sample_id = sample_id,
            merged_vcf_gz = merged_vcf_gz,
            merged_tbi = merged_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai
    }
    
    output {
        File regenotyped_kanpig = GetRegenotypedVcfImpl.regenotyped_kanpig
        File regenotyped_kanpig_tbi = GetRegenotypedVcfImpl.regenotyped_kanpig_tbi
    }
}


task GetRegenotypedVcfImpl {
    input {
        String sample_id
        File merged_vcf_gz
        File merged_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    Int mem_gb = 32
    Int mem_gb_sort = 10
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        KANPIG_SIZEMAX="10000"  # From Adam's suggestion:
        # "I noticed that you used for kanpig --sizemax 1000000 . You're going to get lower recall with that. Currently kanpig is using a very naive clustering strategy to figure out which variants should be considered together. The boundaries of the variant graphs are set to min_start/max_end and pileups are made over the region. Since kanpig is also only looking at pileups of reads that span the region, large variants can preclude smaller variants from getting a chance to have read support. I'm working on better clustering and not needing only spanning reads, but for now sizemax should be set to something like 10k, or maybe even ~75% of the mean insert size."
        chmod +x ~{docker_dir}/kanpig

        # Making sure the VCF has the right sample ID
        bcftools view --header-only ~{merged_vcf_gz} > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > cleaned.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t~{sample_id}" >> cleaned.vcf
        bcftools view --no-header ~{merged_vcf_gz} >> cleaned.vcf
        rm -f ~{merged_vcf_gz}
        bgzip cleaned.vcf
        tabix -f cleaned.vcf.gz

        # Re-genotyping
        rm -f ~{alignments_bai}
        samtools index -@ ${N_THREADS} ~{alignments_bam}
        export RUST_BACKTRACE=1
        ~{docker_dir}/kanpig --threads ${N_THREADS} --sizemin 0 --sizemax ${KANPIG_SIZEMAX} --input cleaned.vcf.gz --bam ~{alignments_bam} --reference ~{reference_fa} --out tmp1.vcf.gz
        bcftools sort --max-mem ~{mem_gb_sort} --output-type z tmp1.vcf.gz > regenotyped_kanpig.vcf.gz
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
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}
