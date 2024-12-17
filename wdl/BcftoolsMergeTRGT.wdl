version 1.0


# Performs a trivial bcftools merge of multiple samples called with TRGT.
# The output VCF contains no multiallelic records and no exact duplicates.
#
workflow BcftoolsMergeTRGT {
    input {
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        Int compression_level = 1
        Int ram_gb = 64
        Int n_cpu = 4
    }
    parameter_meta {
    }
    
    call BcftoolsMergeTRGTImpl {
        input:
            sample_vcf_gz = sample_vcf_gz,
            sample_tbi = sample_tbi,
            compression_level = compression_level,
            ram_gb = ram_gb,
            n_cpu = n_cpu
    }
    
    output {
        File output_vcf_gz = BcftoolsMergeTRGTImpl.output_vcf_gz
        File output_tbi = BcftoolsMergeTRGTImpl.output_tbi
    }
}


task BcftoolsMergeTRGTImpl {
    input {
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
        Int compression_level
        Int ram_gb
        Int n_cpu
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(sample_vcf_gz, "GB")) + 10
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
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        
        function cleanVCF() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF=$2
            
            # Ensuring that the input file is sorted
            ${TIME_COMMAND} bcftools sort --max-mem ${EFFECTIVE_RAM_GB}G --output-type z ${INPUT_VCF_GZ} > tmp0.vcf.gz
            tabix -f tmp0.vcf.gz
            
            # Removing multiallelic records.
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp0.vcf.gz > tmp1.vcf.gz
            tabix -f tmp1.vcf.gz
            rm -f tmp0.vcf.gz*
            
            # Removing identical records
            # See <https://github.com/samtools/bcftools/issues/1089>.
            ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --rm-dup exact --output-type v tmp1.vcf.gz > ${OUTPUT_VCF}
            rm -f tmp1.vcf.gz*
            
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} --compress-level ~{compression_level} ${OUTPUT_VCF}
            tabix -f ${OUTPUT_VCF}.gz
        }
        
        # Merging all single-sample VCFs
        INPUT_FILES=~{sep=',' sample_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        i="0"; SAMPLE_ID="";
        for INPUT_FILE in ${INPUT_FILES}; do
            SAMPLE_ID=$(bcftools view --header-only ${INPUT_FILE} | tail -n 1 | cut -f 10)
            cleanVCF ${INPUT_FILE} ${i}.vcf
            echo ${i}.vcf.gz >> list.txt
            rm -f ${INPUT_FILE}
            i=$(( ${i} + 1 ))
        done
        # $--info-rules -$ disables default rules, and it is used just to avoid
        # the following error:
        # Only fixed-length vectors are supported with -i sum:AC
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --force-samples --file-list list.txt --output-type z > tmp1.vcf.gz
        tabix -f tmp1.vcf.gz
        
        # Removing multiallelic records, if any are generated during the merge.
        # This is just an extra safeguard and might be dropped if $bcftools
        # merge --merge none$ always behaves correctly.
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z tmp1.vcf.gz > out.vcf.gz
        tabix -f out.vcf.gz
        ls -laht; tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: n_cpu
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
