version 1.0


# 
#
workflow BcftoolsMergeIntrasampleBND {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    call BcftoolsMergeIntrasampleBNDImpl {
        input:
            sample_id = sample_id,
            pbsv_vcf_gz = pbsv_vcf_gz,
            pbsv_vcf_gz_tbi = pbsv_vcf_gz_tbi,
            sniffles_vcf_gz = sniffles_vcf_gz,
            sniffles_vcf_gz_tbi = sniffles_vcf_gz_tbi,
            reference_fa = reference_fa
    }
    
    output {
    	File bcftools_merged = BcftoolsMergeIntrasampleBNDImpl.bcftools_merged
    	File bcftools_merged_idx = BcftoolsMergeIntrasampleBNDImpl.bcftools_merged_idx
    }
}


#
task BcftoolsMergeIntrasampleBNDImpl {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(pbsv_vcf_gz,"GB")) + ceil(size(sniffles_vcf_gz,"GB")) ) + 50
    Int mem_gb = 16
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{mem_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        
        # Removing multiallelic records from the input, if any.
        bcftools norm --multiallelics - --output-type z ~{pbsv_vcf_gz} > pbsv_1.vcf.gz
        tabix -f pbsv_1.vcf.gz
        bcftools norm --multiallelics - --output-type z ~{sniffles_vcf_gz} > sniffles_1.vcf.gz
        tabix -f sniffles_1.vcf.gz
        
        # Keeping only BNDs
        bcftools filter --include 'SVTYPE="BND"' --output-type z pbsv_1.vcf.gz > pbsv_2.vcf.gz
        tabix -f pbsv_2.vcf.gz
        bcftools filter --include 'SVTYPE="BND"' --output-type z sniffles_1.vcf.gz > sniffles_2.vcf.gz
        tabix -f sniffles_2.vcf.gz
        rm -f *_1.vcf.gz*
        
        # Fixing REF=N (caused e.g. by sniffles).
        bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z pbsv_2.vcf.gz > pbsv_3.vcf.gz
        tabix -f pbsv_3.vcf.gz
        bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z sniffles_2.vcf.gz > sniffles_3.vcf.gz
        tabix -f sniffles_3.vcf.gz
        rm -f *_2.vcf.gz*

        # Removing exact duplicates
        bcftools concat --threads ${N_THREADS} --allow-overlaps --remove-duplicates --output-type z --output tmp.vcf.gz pbsv_3.vcf.gz sniffles_3.vcf.gz
        tabix -f tmp.vcf.gz
        
        # Removing multiallelic records again, just to be completely sure they
        # don't exist.
        bcftools norm --multiallelics - --output-type z tmp.vcf.gz > ~{sample_id}.bcftools_merged.vcf.gz
        tabix -f ~{sample_id}.bcftools_merged.vcf.gz
        rm -f tmp.vcf.gz*
    >>>
    
    output {
    	File bcftools_merged = "~{work_dir}/~{sample_id}.bcftools_merged.vcf.gz"
    	File bcftools_merged_idx = "~{work_dir}/~{sample_id}.bcftools_merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 2
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
