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
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        Int min_sv_length = 10000
        Int single_breakend_length = 1000
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
            pav_vcf_gz = pav_vcf_gz,
            pav_vcf_gz_tbi = pav_vcf_gz_tbi,
            min_sv_length = min_sv_length,
            single_breakend_length = single_breakend_length
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
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        String remote_chromosomes_dir
        Int min_sv_length
        Int single_breakend_length
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(pbsv_vcf_gz,"GB")) + ceil(size(sniffles_vcf_gz,"GB")) + ceil(size(pav_vcf_gz,"GB")) ) + 50
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
        
        # Remark: the following does not work for BNDs, e.g. it changes ALTs as
        # follows: ]chr2:193700779]N  ->  ccNNNNNNNNNNNNNNN
        ## Fixing REF=N (caused e.g. by sniffles).
        #bcftools norm --check-ref s --fasta-ref reference_fa --do-not-normalize --output-type z pbsv_2.vcf.gz > pbsv_3.vcf.gz
        
        gsutil -m cp ~{remote_chromosomes_dir}/'*' .
        
        # Removing multiallelic records from the input, if any.
        bcftools norm --multiallelics - --output-type z ~{pbsv_vcf_gz} > pbsv_1.vcf.gz
        tabix -f pbsv_1.vcf.gz
        bcftools norm --multiallelics - --output-type z ~{sniffles_vcf_gz} > sniffles_1.vcf.gz
        tabix -f sniffles_1.vcf.gz
        bcftools norm --multiallelics - --output-type v ~{pav_vcf_gz} > pav_1.vcf
        
        # Harvesting the original BNDs from sniffles and pbsv
        bcftools filter --include 'SVTYPE="BND"' --output-type v pbsv_1.vcf.gz > pbsv_2.vcf
        java -cp ~{docker_dir} CleanBNDs pbsv_2.vcf . pbsv_bnds.vcf
        bcftools sort --output-type z pbsv_bnds.vcf > pbsv_bnds.vcf.gz
        tabix -f pbsv_bnds.vcf.gz
        bcftools filter --include 'SVTYPE="BND"' --output-type v sniffles_1.vcf.gz > sniffles_2.vcf
        java -cp ~{docker_dir} CleanBNDs sniffles_2.vcf . sniffles_bnds.vcf
        bcftools sort --output-type z sniffles_bnds.vcf > sniffles_bnds.vcf.gz
        tabix -f sniffles_bnds.vcf.gz
        rm -f *_2.vcf
        
        # Creating new BNDs from the large calls of every caller
        bcftools filter --include 'SVTYPE!="BND"' --output-type v pbsv_1.vcf.gz > pbsv_2.vcf
        java -cp ~{docker_dir} SV2BND pbsv_2.vcf . ~{min_sv_length} ~{single_breakend_length} pbsv_others.vcf
        bcftools sort --output-type z pbsv_others.vcf > pbsv_others.vcf.gz
        tabix -f pbsv_others.vcf.gz
        bcftools filter --include 'SVTYPE!="BND"' --output-type v sniffles_1.vcf.gz > sniffles_2.vcf
        java -cp ~{docker_dir} SV2BND sniffles_2.vcf . ~{min_sv_length} ~{single_breakend_length} sniffles_others.vcf
        bcftools sort --output-type z sniffles_others.vcf > sniffles_others.vcf.gz
        tabix -f sniffles_others.vcf.gz
        java -cp ~{docker_dir} SV2BND pav_1.vcf . ~{min_sv_length} ~{single_breakend_length} pav_others.vcf
        bcftools sort --output-type z pav_others.vcf > pav_others.vcf.gz
        tabix -f pav_others.vcf.gz
        rm -f *_2.vcf
        
        # Removing exact duplicates
        N_BNDS_PBSV=$(bcftools view --no-header pbsv_bnds.vcf.gz | wc -l)
        N_BNDS_SNIFFLES=$(bcftools view --no-header sniffles_bnds.vcf.gz | wc -l)
        N_OTHERS_PBSV=$(bcftools view --no-header pbsv_others.vcf.gz | wc -l)
        N_OTHERS_SNIFFLES=$(bcftools view --no-header sniffles_others.vcf.gz | wc -l)
        N_OTHERS_PAV=$(bcftools view --no-header pav_others.vcf.gz | wc -l)
        bcftools concat --threads ${N_THREADS} --allow-overlaps --remove-duplicates --output-type z --output tmp.vcf.gz pbsv_bnds.vcf.gz sniffles_bnds.vcf.gz pbsv_others.vcf.gz sniffles_others.vcf.gz pav_others.vcf.gz
        tabix -f tmp.vcf.gz
        
        # Removing multiallelic records again, just to be completely sure they
        # don't exist.
        bcftools norm --multiallelics - --output-type z tmp.vcf.gz > ~{sample_id}.bcftools_merged.vcf.gz
        tabix -f ~{sample_id}.bcftools_merged.vcf.gz
        rm -f tmp.vcf.gz*
        N_BNDS_FINAL=$(bcftools view --no-header ~{sample_id}.bcftools_merged.vcf.gz | wc -l)
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
