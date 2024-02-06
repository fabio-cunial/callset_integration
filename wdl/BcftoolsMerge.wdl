version 1.0


#
workflow BcftoolsMerge {
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
        call IntraSampleConcat {
            input:
                sample_vcf_gz = [pbsv_vcf_gz[i], sniffles_vcf_gz[i], pav_vcf_gz[i]],
                sample_tbi = [pbsv_tbi[i], sniffles_tbi[i], pav_tbi[i]]
        }
    }
    call InterSampleMerge {
        input:
            input_vcf_gz = IntraSampleConcat.output_vcf_gz,
            input_tbi = IntraSampleConcat.output_tbi
    }
    
    output {
        File output_vcf_gz = InterSampleMerge.output_vcf_gz
        File output_tbi = InterSampleMerge.output_tbi
    }
}


# Remark: $bcftools concat$ does not try to merge any records, so the output
# of this task contains all the original records.
#
task IntraSampleConcat {
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
        
        INPUT_FILES=~{sep=',' sample_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list.txt --output-type z > concat.vcf.gz
        tabix concat.vcf.gz
        ls -laht; tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/concat.vcf.gz"
        File output_tbi = work_dir + "/concat.vcf.gz.tbi"
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
# differences in SVLEN, END or STRAND. This may delete information for symbolic
# ALTs. Our script makes sure that only symbolic records with the same SVLEN,
# END and STRAND are collapsed into the same record.
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
                    strand=""; \
                    n=split($8,A,";"); \
                    for (i=1; i<=n; i++) { \
                        if (substr(A[i],1,4)=="END=") end=substr(A[i],5); \
                        else if (substr(A[i],1,6)=="SVLEN=") { \
                            if (substr(A[i],7,1)=="-") svlen=substr(A[i],8); \
                            else svlen=substr(A[i],7); \
                        } \
                        else if (substr(A[i],1,7)=="STRAND=") strand=substr(A[i],8); \
                    } \
                    $5="<" svtype ":" tag ":" (length(end)==0?"?":end) ":" (length(svlen)==0?"?":svlen) ":" (length(strand)==0?"?":strand) ">" \
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
            ID="$(basename ${INPUT_FILE} .vcf.gz)_${i}"
            cleanVCF ${INPUT_FILE} ${ID}_cleaned.vcf
            echo ${ID}_cleaned.vcf.gz >> list.txt
            rm -f ${INPUT_FILE}
            i=$(( ${i} + 1 ))
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --file-list list.txt --output-type v > tmp.vcf
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
