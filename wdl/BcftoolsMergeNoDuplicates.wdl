version 1.0


#
workflow BcftoolsMergeNoDuplicates {
    input {
        Array[File] pbsv_vcf_gz
        Array[File] pbsv_tbi
        Array[File] sniffles_vcf_gz
        Array[File] sniffles_tbi
        Array[File] pav_vcf_gz
        Array[File] pav_tbi
    }
    parameter_meta {
    }
    
    scatter(i in range(length(pbsv_vcf_gz))) {
        call IntraSampleMerge {
            input:
                sample_vcf_gz = [pbsv_vcf_gz[i], sniffles_vcf_gz[i], pav_vcf_gz[i]],
                sample_tbi = [pbsv_tbi[i], sniffles_tbi[i], pav_tbi[i]]
        }
    }
    call InterSampleMerge {
        input:
            input_vcf_gz = IntraSampleMerge.output_vcf_gz,
            input_tbi = IntraSampleMerge.output_tbi
    }
    
    output {
        File output_vcf_gz = InterSampleMerge.output_vcf_gz
        File output_tbi = InterSampleMerge.output_tbi
    }
}


# Remark: we use $bcftools merge$ instead of $bcftools concat$, since we must
# collapse identical calls made by different callers (otherwise they would
# remain in the inter-sample VCF, since the inter-sample $bcftools merge$ does
# not collapse records from the same sample).
#
# Remark: $bcftools merge$ collapses into the same record every record with the
# same $CHR,POS,REF,ALT$, disregarding the INFO field and in particular 
# differences in SVLEN and END. This may delete information for symbolic
# ALTs. Our script makes sure that only symbolic records with the same SVLEN and
# END are collapsed into the same record.
#
# Remark: we do not consider STRAND in the above.
#
# Remark: the output VCF has artificial FORMAT and SAMPLE columns where all
# calls are 0/1. I.e. the original genotypes are discarded.
#
task IntraSampleMerge {
    input {
        Array[File] sample_vcf_gz
        Array[File] sample_tbi
    }
    parameter_meta {
    }
    
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
        
        function cleanVCF() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF=$2
            
            bcftools view --header-only ${INPUT_VCF_GZ} > ${OUTPUT_VCF}
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF_GZ} | awk '{ \
                tag="artificial"; \
                if ($5=="<DEL>" || $5=="<INS>" || $5=="<INV>" || $5=="<DUP>" || $5=="<CNV>") { \
                    svtype=substr($5,2,3); \
                    end=""; \
                    svlen=""; \
                    n=split($8,A,";"); \
                    for (i=1; i<=n; i++) { \
                        if (substr(A[i],1,4)=="END=") end=substr(A[i],5); \
                        else if (substr(A[i],1,6)=="SVLEN=") { \
                            if (substr(A[i],7,1)=="-") svlen=substr(A[i],8); \
                            else svlen=substr(A[i],7); \
                        } \
                    } \
                    $5="<" svtype ":" tag ":" (length(end)==0?"?":end) ":" (length(svlen)==0?"?":svlen) ">" \
                }; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> ${OUTPUT_VCF}
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} ${OUTPUT_VCF}
            tabix ${OUTPUT_VCF}.gz
        }
        
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
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > tmp1.vcf.gz
        bcftools view --header-only tmp1.vcf.gz > header.txt
        N_ROWS=$(wc -l < header.txt)
        head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
        echo -e "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t${SAMPLE_ID}" >> out.vcf
        ${TIME_COMMAND} bcftools view --no-header tmp1.vcf.gz | awk '{ \
            tag="artificial"; \
            if (substr($5,6,length(tag))==tag) $5=substr($5,1,4) ">"; \
            printf("%s",$1); \
            for (i=2; i<=8; i++) printf("\t%s",$i); \
            printf("\tGT\t0/1\n"); \
        }' >> out.vcf
        bgzip out.vcf
        tabix -f out.vcf.gz
        ls -laht; tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/out.vcf.gz"
        File output_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: "8GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}


# Remark: $bcftools merge$ collapses into the same record every record with the
# same $CHR,POS,REF,ALT$, disregarding the INFO field and in particular 
# differences in SVLEN and END. This may delete information for symbolic
# ALTs. Our script makes sure that only symbolic records with the same SVLEN and
# END are collapsed into the same record.
#
# Remark: we do not consider STRAND in the above.
#
# Remark: since the input comes from $IntraSampleMerge$, which overwrites GTs,
# every call in the output of this procedure has at least one sample with
# GT=0/1, and every sample has GT \in {0/1, ./.}.
#
task InterSampleMerge {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
    }
    parameter_meta {
    }
    
    Int disk_size_gb = ceil(size(input_vcf_gz, "GB")) + 100
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
        
        function cleanVCF() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF=$2
            
            bcftools view --header-only ${INPUT_VCF_GZ} > ${OUTPUT_VCF}
            ${TIME_COMMAND} bcftools view --no-header ${INPUT_VCF_GZ} | awk '{ \
                tag="artificial"; \
                if ($5=="<DEL>" || $5=="<INS>" || $5=="<INV>" || $5=="<DUP>" || $5=="<CNV>") { \
                    svtype=substr($5,2,3); \
                    end=""; \
                    svlen=""; \
                    n=split($8,A,";"); \
                    for (i=1; i<=n; i++) { \
                        if (substr(A[i],1,4)=="END=") end=substr(A[i],5); \
                        else if (substr(A[i],1,6)=="SVLEN=") { \
                            if (substr(A[i],7,1)=="-") svlen=substr(A[i],8); \
                            else svlen=substr(A[i],7); \
                        } \
                    } \
                    $5="<" svtype ":" tag ":" (length(end)==0?"?":end) ":" (length(svlen)==0?"?":svlen) ">" \
                }; \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            }' >> ${OUTPUT_VCF}
            ${TIME_COMMAND} bgzip --threads ${N_THREADS} ${OUTPUT_VCF}
            tabix ${OUTPUT_VCF}.gz
        }
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        i="0"
        for INPUT_FILE in ${INPUT_FILES}; do
            cleanVCF ${INPUT_FILE} ${i}.vcf
            echo ${i}.vcf.gz >> list.txt
            rm -f ${INPUT_FILE}
            i=$(( ${i} + 1 ))
        done
        # $--info-rules -$ disables default rules, and it is used just to avoid
        # the following error: 
        # Only fixed-length vectors are supported with -i sum:AC
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --info-rules - --file-list list.txt --output-type v > tmp.vcf
        bcftools view --header-only tmp.vcf > merged.vcf
        ${TIME_COMMAND} bcftools view --no-header tmp.vcf | awk '{ \
            tag="artificial"; \
            if (substr($5,6,length(tag))==tag) $5=substr($5,1,4) ">"; \
            printf("%s",$1); \
            for (i=2; i<=NF; i++) printf("\t%s",$i); \
            printf("\n"); \
        }' >> merged.vcf
        rm -f tmp.vcf
        ${TIME_COMMAND} bgzip --threads ${N_THREADS} merged.vcf
        tabix merged.vcf.gz
        ls -laht; tree
    >>>
    
    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: "32GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
