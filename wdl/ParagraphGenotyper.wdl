version 1.0


# Features added by Paragraph: GT:GQ:DR:DV:OLD_GT:DP:FT:AD:ADF:ADR:PL
#
# INFO.GRMPY_ID="Graph ID for linking to genotypes.json.gz; matches
# record.graphinfo.ID in there.
#
# FORMAT.OLD_GT="Previous GT which was replaced by paragraph"
# FORMAT.FT="Filter for genotype"
# FORMAT.ADF="Allele depth on forward strand for each allele, including the
# reference."
# FORMAT.ADR="Allele depth on reverse strand for each allele, including the
# reference."
# FORMAT.PL="Phred-scaled likelihoods for genotypes as defined in the VCF
# specification"
#
# FILTER.BP_DEPTH="One or more breakpoints have abnormal depth"
# FILTER.NO_VALID_GT="No valid genotypes from breakpoints"
# FILTER.CONFLICT="Breakpoints gave different genotypes"
# FILTER.BP_NO_GT="One genotype was missing"
# FILTER.NO_READS="No reads could be retrieved for a breakpoint."
# FILTER.DEPTH="Poisson depth filter: observed depth deviates too far from
# Poisson expectation"
# FILTER.UNMATCHED="VCF record could not be matched to a paragraph record."
# FILTER.MULTIMATCHED="VCF record could not be matched to a paragraph record
# uniquely."
#
workflow ParagraphGenotyper {
    input {
        File input_vcf_gz
        File short_reads_bam
        File short_reads_bai
        Int short_read_coverage
        Int short_read_length
        File reference_fa
        File reference_fai
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    call ParagraphImpl {
        input:
        input_vcf_gz = input_vcf_gz,
        short_reads_bam = short_reads_bam,
        short_reads_bai = short_reads_bai,
        short_read_coverage = short_read_coverage,
        short_read_length = short_read_length,
        reference_fa = reference_fa,
        reference_fai = reference_fai,
        n_cpus = n_cpus,
        ram_size_gb = ram_size_gb,
        disk_size_gb = disk_size_gb
    }
    
    output {
        File out_vcf_gz = ParagraphImpl.out_vcf_gz
        File out_tbi = ParagraphImpl.out_tbi
    }
}


# TIME | CORES | RAM   | COST
# 1.5h |  27   | 6 GB  | 2.50 $
#
task ParagraphImpl {
    input {
        File input_vcf_gz
        File short_reads_bam
        File short_reads_bai
        Int short_read_coverage
        Int short_read_length
        File reference_fa
        File reference_fai
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        PARAGRAPH_DIR="~{docker_dir}/paragraph/bin"
        
        # Cleaning paragraph's input
        gunzip -c ~{input_vcf_gz} > input.vcf
        rm -f ~{input_vcf_gz}
        java -cp ~{docker_dir} ParagraphClean input.vcf clean.vcf
        rm -f input.vcf
        
        # Running paragraph
        echo -e "id\tpath\tdepth\tread length" > manifest.txt
        echo -e "sample1\t~{short_reads_bam}\t~{short_read_coverage}\t~{short_read_length}" >> manifest.txt
        mkdir ./tmpdir
        cd ${PARAGRAPH_DIR}
        ls -laht
        ${TIME_COMMAND} python3 multigrmpy.py --threads ${N_THREADS} --verbose --input ~{work_dir}/clean.vcf --reference-sequence ~{reference_fa} --manifest ~{work_dir}/manifest.txt --output ~{work_dir}/tmpdir
        ls -laht
        
        # Formatting paragraph's output
        cd ~{work_dir}
        bcftools view -h ./tmpdir/genotypes.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> out.vcf
        bcftools view -H ./tmpdir/genotypes.vcf.gz | cut -f 1-9,11 >> out.vcf
        bgzip out.vcf
        tabix -f out.vcf.gz
    >>>
    
    output {
        File out_vcf_gz = work_dir + "/out.vcf.gz"
        File out_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
