version 1.0


# 
#
workflow BcftoolsMergeIntrasample {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    call BcftoolsMergeIntrasampleImpl {
        input:
            sample_id = sample_id,
            pbsv_vcf_gz = pbsv_vcf_gz,
            pbsv_vcf_gz_tbi = pbsv_vcf_gz_tbi,
            sniffles_vcf_gz = sniffles_vcf_gz,
            sniffles_vcf_gz_tbi = sniffles_vcf_gz_tbi,
            pav_vcf_gz = pav_vcf_gz,
            pav_vcf_gz_tbi = pav_vcf_gz_tbi,
            reference_fa = reference_fa
    }
    
    output {
    	File bcftools_merged = BcftoolsMergeIntrasampleImpl.bcftools_merged
    	File bcftools_merged_idx = BcftoolsMergeIntrasampleImpl.bcftools_merged_idx
    }
}


#
task BcftoolsMergeIntrasampleImpl {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*( ceil(size(pbsv_vcf_gz,"GB")) + ceil(size(sniffles_vcf_gz,"GB")) + ceil(size(pav_vcf_gz,"GB")) + ceil(size(reference_fa,"GB")) ) + 50
    Int mem_gb = 128
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
        tabix pbsv_1.vcf.gz
        bcftools norm --multiallelics - --output-type z ~{sniffles_vcf_gz} > sniffles_1.vcf.gz
        tabix sniffles_1.vcf.gz
        bcftools norm --multiallelics - --output-type z ~{pav_vcf_gz} > pav_1.vcf.gz
        tabix pav_1.vcf.gz
        
        # Fixing REF=N (caused e.g. by sniffles).
        bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z pbsv_1.vcf.gz > pbsv_new.vcf.gz
        tabix -f pbsv_new.vcf.gz
        bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z sniffles_1.vcf.gz > sniffles_new.vcf.gz
        tabix -f sniffles_new.vcf.gz
        bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z pav_1.vcf.gz > pav_new.vcf.gz
        tabix -f pav_new.vcf.gz
        rm -f *_1.vcf.gz*
        
        # Putting every file in a consistent format
        mkdir -p preprocessed
        for in_vcf in pav_new.vcf.gz pbsv_new.vcf.gz sniffles_new.vcf.gz
        do
            outname=preprocessed/$(basename $in_vcf)
            python ~{docker_dir}/resolve_light.py ${in_vcf} ~{reference_fa} | bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z > ${outname}
            tabix -f $outname
        done

        # Removing exact duplicates
        bcftools concat --threads ${N_THREADS} --allow-overlaps --remove-duplicates --output-type z --output tmp.vcf.gz preprocessed/pbsv_new.vcf.gz preprocessed/sniffles_new.vcf.gz preprocessed/pav_new.vcf.gz
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
