version 1.0


# 
#
workflow MergeSVsSNPs {
    input {
        File svs_vcf_gz
        File svs_tbi
        File snps_vcf_gz
        File snps_tbi
        File reference_fa
        File reference_fai
        Int ram_size_gb = 256
    }
    parameter_meta {
    }
    
    call MergeSVsSNPsImpl {
        input:
            svs_vcf_gz = svs_vcf_gz,
            svs_tbi = svs_tbi,
            snps_vcf_gz = snps_vcf_gz,
            snps_tbi = snps_tbi,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            ram_size_gb = ram_size_gb
    }
    
    output {
        File output_vcf_gz = MergeSVsSNPsImpl.output_vcf_gz
        File output_tbi = MergeSVsSNPsImpl.output_tbi
    }
}


task MergeSVsSNPsImpl {
    input {
        File svs_vcf_gz
        File svs_tbi
        File snps_vcf_gz
        File snps_tbi
        File reference_fa
        File reference_fai
        Int ram_size_gb
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(svs_vcf_gz, "GB")+size(snps_vcf_gz, "GB")) + 10
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_size_gb} - 2 ))
        
        # Cleaning the SNP VCF:
        # - removing multiallelic records;
        # - fixing wrong REF values.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z ~{snps_vcf_gz} > snps.vcf.gz
        rm -f ~{snps_vcf_gz}
        ls -laht; tree; df -h
        tabix -f snps.vcf.gz
        
        # Merging SNPs and SVs
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --remove-duplicates --output-type z ~{svs_vcf_gz} snps.vcf.gz > merged.vcf.gz
        ls -laht; tree; df -h
        tabix -f merged.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
