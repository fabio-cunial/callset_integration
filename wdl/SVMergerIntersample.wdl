version 1.0


# 
#
workflow SVMergerIntersample {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        String remote_trf_dir
    }
    parameter_meta {
        input_vcf_gz: "One VCF per sample (from intra-sample merging)."
        remote_trf_dir: "Contains per-chromosome TRF files."
    }
    
    call GetChromosomeTSVs {
        input:
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi
    }
    Array[File] chromosome_vcf_gz = GetChromosomeTSVs.chromosome_vcf_gz
    Array[File] chromosome_tbi = GetChromosomeTSVs.chromosome_tbi
    Array[File] chromosome_tsv = GetChromosomeTSVs.chromosome_tsv
    scatter(i in range(length(chromosome_tsv))) {
        call SVMergerImpl {
            input:
                input_vcf_gz = chromosome_vcf_gz[i],
                input_tbi = chromosome_tbi[i],
                input_tsv = chromosome_tsv[i],
                remote_trf_dir = remote_trf_dir
        }
    }
    call ConcatChromosomeVCFs {
        input:
            chromosome_vcf_gz = SVMergerImpl.merged_vcf_gz,
            chromosome_tbi = SVMergerImpl.merged_tbi
    }
    
    output {
        File output_vcf_gz = ConcatChromosomeVCFs.output_vcf_gz
        File output_vcf_gz_tbi = ConcatChromosomeVCFs.output_tbi
    }
}


# Given all single-sample VCFs from all samples, the task outputs one
# concatenated TSV and one bcftools-merged VCF per chromosome.
#
task GetChromosomeTSVs {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
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
        
        cp ~{docker_dir}/*.class .
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        
        # 1. Merging all VCFs and splitting the merge by chromosome.
        # 1.1 Calls from different samples might have the same ID, which must be
        # made unique.
        # 1.2 $bcftools merge$ might merge identical calls from different 
        # samples: their IDs are concatenated in the merged file.
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            SAMPLE_ID=$(basename ${INPUT_FILE} .vcf.gz)
            ${TIME_COMMAND} bcftools annotate --threads ${N_THREADS} --set-id ${SAMPLE_ID}'-%ID' ${INPUT_FILE} --output-type z > ${SAMPLE_ID}-annotated.vcf.gz
            tabix ${SAMPLE_ID}-annotated.vcf.gz
            echo ${SAMPLE_ID}-annotated.vcf.gz >> list.txt
            rm -f ${INPUT_FILE}
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > merge.vcf.gz
        tabix merge.vcf.gz
        rm -f output_vcf_gz_list.txt output_tbi_list.txt
        for CHR in $(seq 1 22) X Y; do
            bcftools view --output-type z merge.vcf.gz chr${CHR} > chr${CHR}.vcf.gz
            tabix chr${CHR}.vcf.gz
            echo ~{work_dir}/chr${CHR}.vcf.gz >> output_vcf_gz_list.txt
            echo ~{work_dir}/chr${CHR}.vcf.gz.tbi >> output_tbi_list.txt
        done
        rm -f merge.vcf.gz
        
        # 2. Creating per-chromosome TSVs.
        # Remark: we must run $VCF2SVMerger$ on the output of bcftools merge, to
        # make sure the TSV uses the same IDs as the VCF.
        rm -f output_tsv_list.txt
        for CHR in $(seq 1 22) X Y; do
            gunzip --stdout chr${CHR}.vcf.gz > tmp.vcf
            ${TIME_COMMAND} java VCF2SVMerger tmp.vcf ${SAMPLE_ID} null > chr${CHR}.tsv
            echo ~{work_dir}/chr${CHR}.tsv >> output_tsv_list.txt
            rm -f tmp.vcf
        done
        ls -laht; tree
    >>>
    
    output {
        # Using glob() instead of read_line() since reading files from a list is
        # not supported. I hope different calls to glob() return files in the
        # same order.
        # https://support.terra.bio/hc/en-us/community/posts/21977172563099-Specify-Array-File-output-for-a-task-in-a-WDL
        Array[File] chromosome_vcf_gz = glob(work_dir+"/chr*.vcf.gz")
        Array[File] chromosome_tbi = glob(work_dir+"/chr*.vcf.gz.tbi")
        Array[File] chromosome_tsv = glob(work_dir+"/chr*.tsv")
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: "16GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}


# Given a chromosome TSV and a chromosome merged VCF, the task runs svmerger
# and outputs its result as a VCF.
#
task SVMergerImpl {
    input {
        File input_vcf_gz
        File input_tbi
        File input_tsv
        String remote_trf_dir
    }
    parameter_meta {
        input_vcf_gz: "All calls from all samples for a given chromosome."
        input_tsv: "All calls from all samples for a given chromosome, encoded in TSV format."
        remote_trf_dir: "Contains per-chromosome TRF files."
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
        
        # 1. Running svmerger on the input TSV to get clique representatives
        cp ~{docker_dir}/*.class .
        CHROMOSOME_ID=$(basename ~{input_tsv} .tsv)
        gsutil -m cp ~{remote_trf_dir}/${CHROMOSOME_ID}.trf.sorted.gor  .
        
        function svMerger() {
            local CHROMOSOME_ID=$1
            local SVTYPE=$2
            
            if [ ${SVTYPE} = "DEL" ]; then
                AWK_COMMAND='$7 == "DEL"'
            elif [ ${SVTYPE} = "INS" ]; then
                AWK_COMMAND='$7 == "INS"'
            elif [ ${SVTYPE} = "DUP" ]; then
                AWK_COMMAND='$7 == "DUP"'
            elif [ ${SVTYPE} = "INV" ]; then
                AWK_COMMAND='$7 == "INV"'
            fi
            awk "${AWK_COMMAND}" ~{input_tsv} > ${CHROMOSOME_ID}.${SVTYPE}.tsv
            if [ -s ${CHROMOSOME_ID}.${SVTYPE}.tsv ]; then
                ${TIME_COMMAND} python2 ~{docker_dir}/sv-merger/main.py MERGE ${CHROMOSOME_ID}.${SVTYPE}.tsv ${CHROMOSOME_ID}.trf.sorted.gor ${SVTYPE}
                sort -k 4 ${CHROMOSOME_ID}.${SVTYPE}.tsv > chr.tsv
                ls -laht; tree
                # Outside TRs
                if [ -s ${CHROMOSOME_ID}.${SVTYPE}.tsv.outtrr.merged.csv ]; then
                    sort -k 1 ${CHROMOSOME_ID}.${SVTYPE}.tsv.outtrr.merged.csv > outtrr.tsv
                    join -t $'\t' -1 4 -2 1 chr.tsv outtrr.tsv | sort --version-sort --key 9 > cliques.tsv
                    ls -laht; tree
                    java SVMergerGetRepresentative cliques.tsv > clique-representatives-${SVTYPE}-outttr.tsv
                    rm -f cliques.tsv outtrr.tsv ${CHROMOSOME_ID}.${SVTYPE}.tsv.outtrr.merged.csv
                fi
                # Inside TRs
                if [ -s ${CHROMOSOME_ID}.${SVTYPE}.tsv.intrr.merged.csv ]; then
                    sort -k 1 ${CHROMOSOME_ID}.${SVTYPE}.tsv.intrr.merged.csv > intrr.tsv
                    join -t $'\t' -1 4 -2 1 chr.tsv intrr.tsv | sort --version-sort --key 9 > cliques.tsv
                    ls -laht; tree
                    java SVMergerGetRepresentative cliques.tsv > clique-representatives-${SVTYPE}-inttr.tsv
                    rm -f cliques.tsv intrr.tsv ${CHROMOSOME_ID}.${SVTYPE}.tsv.intrr.merged.csv
                fi
            fi
        }
        
        svMerger ${CHROMOSOME_ID} DEL
        ls -laht; tree
        svMerger ${CHROMOSOME_ID} INS
        ls -laht; tree
        svMerger ${CHROMOSOME_ID} DUP
        ls -laht; tree
        svMerger ${CHROMOSOME_ID} INV
        ls -laht; tree
        cat clique-representatives-*.tsv | sort > clique-representatives.tsv
        rm -f clique-representatives-*.tsv
        
        # 2. Building a VCF file that contains only clique representatives
        bcftools view --no-header ~{input_vcf_gz} | sort -t$'\t' -k3,3 > input.vcf
        bcftools view --header-only ~{input_vcf_gz} > output.vcf
        join -t$'\t' -1 3 -2 1 input.vcf clique-representatives.tsv | awk 'BEGIN {FS="\t"; OFS="\t"} { printf("%s\t%s\t%s\t",$2,$3,$1); for (i=4; i<NF; i++) printf("%s\t",$i); print $NF }' >> output.vcf
        bcftools sort output.vcf --output-type z > output.vcf.gz
        tabix output.vcf.gz
        rm -f output.vcf
    >>>
    
    output {
        File merged_vcf_gz = work_dir+"/output.vcf.gz"
        File merged_tbi = work_dir+"/output.vcf.gz.tbi"
        File clique_representatives = work_dir+"/clique-representatives.tsv"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: "16GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}


# Creates a single merged VCF from the concatenation of all single-chromosome
# VCFs produced by svmerger.
#
task ConcatChromosomeVCFs {
    input {
        Array[File] chromosome_vcf_gz
        Array[File] chromosome_tbi
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
        
        INPUT_FILES=~{sep=',' chromosome_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list.txt --output-type z > concat.vcf.gz
        tabix concat.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/concat.vcf.gz"
        File output_tbi = work_dir + "/concat.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 4
        memory: "16GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}