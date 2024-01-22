version 1.0


# 
#
workflow SVMerger {
    input {
        String remote_dir
        String sample_id
        String remote_trf_dir
    }
    parameter_meta {
        remote_dir: "Root dir (contains every sample subdir)."
        remote_trf_dir: "Contains per-chromosome TRF files."
    }
    
    call GetTSVs {
        input:
            remote_dir = remote_dir,
            sample_id = sample_id
    }
    Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    scatter(chr in chromosomes) {
        call SVMergerImpl {
            input:
                remote_dir = remote_dir,
                sample_id = sample_id,
                remote_trf_dir = remote_trf_dir,
                chromosome_id = chr,
                artificial_input = GetTSVs.artificial_output
        }
    }
    call SVMerger2VCF {
        input:
            remote_dir = remote_dir,
            sample_id = sample_id,
            artificial_input = SVMergerImpl.artificial_output
    }
    
    output {
        File output_vcf_gz = SVMerger2VCF.output_vcf_gz
        File output_vcf_gz_tbi = SVMerger2VCF.output_vcf_gz_tbi
    }
}


task GetTSVs {
    input {
        String remote_dir
        String sample_id
    }
    parameter_meta {
        remote_dir: "Root dir (contains every sample subdir). Per-chromosome TSVs are stored in the '/tsvs/' subsdirectory of the sample."
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
        rm -f *.tsv
        
        gsutil cp ~{remote_dir}/~{sample_id}/~{sample_id}.pbsv.vcf.gz .
        gunzip --stdout ~{sample_id}.pbsv.vcf.gz > tmp.vcf
        java VCF2SVMerger tmp.vcf ~{sample_id} pbsv >> pbsv.tsv
        rm -f tmp.vcf
        
        gsutil cp ~{remote_dir}/~{sample_id}/~{sample_id}.sniffles.vcf.gz .
        gunzip --stdout ~{sample_id}.sniffles.vcf.gz > tmp.vcf
        java VCF2SVMerger tmp.vcf ~{sample_id} sniffles >> sniffles.tsv
        rm -f tmp.vcf
        
        gsutil cp ~{remote_dir}/~{sample_id}/~{sample_id}.pav_sv.vcf.gz .
        gunzip --stdout ~{sample_id}.pav_sv.vcf.gz > tmp.vcf
        java VCF2SVMerger tmp.vcf ~{sample_id} pav >> pav.tsv
        rm -f tmp.vcf
        
        for CHR in $(seq 1 22) X Y; do
            grep ^chr${CHR}$'\t' pbsv.tsv >> chr${CHR}.tsv
            grep ^chr${CHR}$'\t' sniffles.tsv >> chr${CHR}.tsv
            grep ^chr${CHR}$'\t' pav.tsv >> chr${CHR}.tsv
        done
        gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp '*.tsv' ~{remote_dir}/~{sample_id}/tsvs/
    >>>
    
    output {
        Int artificial_output = 0
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "8GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}


task SVMergerImpl {
    input {
        String remote_dir
        String sample_id
        String remote_trf_dir
        String chromosome_id
        Int artificial_input
    }
    parameter_meta {
        remote_dir: "Root dir (contains every sample subdir)."
        remote_trf_dir: "Contains per-chromosome TRF files."
        artificial_input: "Just to force task execution order."
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
            awk "${AWK_COMMAND}" ${CHROMOSOME_ID}.tsv > ${CHROMOSOME_ID}.${SVTYPE}.tsv
            if [ -s ${CHROMOSOME_ID}.${SVTYPE}.tsv ]; then
                ${TIME_COMMAND} python2 ~{docker_dir}/sv-merger/main.py MERGE ${CHROMOSOME_ID}.${SVTYPE}.tsv ${CHROMOSOME_ID}.trf.sorted.gor ${SVTYPE}
                sort -k 4 ${CHROMOSOME_ID}.${SVTYPE}.tsv > chr.tsv
                ls -laht; tree
                # Outside TRs
                if [ -s ${CHROMOSOME_ID}.${SVTYPE}.tsv.outtrr.merged.csv ]; then
                    sort -k 1 ${CHROMOSOME_ID}.${SVTYPE}.tsv.outtrr.merged.csv > outtrr.tsv
                    join -t $'\t' -1 4 -2 1 chr.tsv outtrr.tsv | sort --version-sort --key 9 > cliques.tsv
                    # tr ' ' '\t' | 
                    ls -laht; tree
                    java SVMergerGetRepresentative cliques.tsv > clique-representatives-${SVTYPE}-outttr.tsv
                    rm -f cliques.tsv outtrr.tsv ${CHROMOSOME_ID}.${SVTYPE}.tsv.outtrr.merged.csv
                fi
                # Inside TRs
                if [ -s ${CHROMOSOME_ID}.${SVTYPE}.tsv.intrr.merged.csv ]; then
                    sort -k 1 ${CHROMOSOME_ID}.${SVTYPE}.tsv.intrr.merged.csv > intrr.tsv
                    join -t $'\t' -1 4 -2 1 chr.tsv intrr.tsv | sort --version-sort --key 9 > cliques.tsv
                    # | tr ' ' '\t'
                    ls -laht; tree
                    java SVMergerGetRepresentative cliques.tsv > clique-representatives-${SVTYPE}-inttr.tsv
                    rm -f cliques.tsv intrr.tsv ${CHROMOSOME_ID}.${SVTYPE}.tsv.intrr.merged.csv
                fi
            fi
        }
        
        cp ~{docker_dir}/*.class .
        gsutil -m cp ~{remote_dir}/~{sample_id}/tsvs/~{chromosome_id}.tsv .
        gsutil -m cp ~{remote_trf_dir}/~{chromosome_id}.trf.sorted.gor  .        
        svMerger ~{chromosome_id} DEL
        ls -laht; tree
        svMerger ~{chromosome_id} INS
        ls -laht; tree
        svMerger ~{chromosome_id} DUP
        ls -laht; tree
        svMerger ~{chromosome_id} INV
        ls -laht; tree
        cat clique-representatives-*.tsv | sort > clique-representatives-~{chromosome_id}.tsv
        gsutil -m cp clique-representatives-~{chromosome_id}.tsv ~{remote_dir}/~{sample_id}/tsvs/
    >>>
    
    output {
        Int artificial_output = 0
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "8GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}


task SVMerger2VCF {
    input {
        String remote_dir
        String sample_id
        Array[Int] artificial_input
    }
    parameter_meta {
        remote_dir: "Root dir (contains every sample subdir)."
        artificial_input: "Just to force task execution order."
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
        
        # Merging all VCFs
        gsutil cp ~{remote_dir}/~{sample_id}/~{sample_id}.pbsv.vcf.'gz*' .
        bcftools annotate --set-id 'pbsv-%ID' ~{sample_id}.pbsv.vcf.gz --output-type z > pbsv.annotated.vcf.gz
        tabix pbsv.annotated.vcf.gz
        rm -f ~{sample_id}.pbsv.vcf.gz*
        gsutil cp ~{remote_dir}/~{sample_id}/~{sample_id}.sniffles.vcf.'gz*' .
        bcftools annotate --set-id 'sniffles-%ID' ~{sample_id}.sniffles.vcf.gz --output-type z > sniffles.annotated.vcf.gz
        tabix sniffles.annotated.vcf.gz
        rm -f ~{sample_id}.sniffles.vcf.gz*
        gsutil cp ~{remote_dir}/~{sample_id}/~{sample_id}.pav_sv.vcf.'gz*' .
        bcftools annotate --set-id 'pav-%ID' ~{sample_id}.pav_sv.vcf.gz --output-type z > pav.annotated.vcf.gz
        tabix pav.annotated.vcf.gz
        rm -f ~{sample_id}.pav_sv.vcf.gz*
        bcftools concat --allow-overlaps --output-type z pbsv.annotated.vcf.gz sniffles.annotated.vcf.gz pav.annotated.vcf.gz > all_calls.vcf.gz
        tabix all_calls.vcf.gz
        rm -f pbsv.annotated.vcf.gz sniffles.annotated.vcf.gz pav.annotated.vcf.gz
        
        # Building a VCF file that contains only clique representatives
        rm -f list.txt
        for CHR in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY ; do
            bcftools view --no-header all_calls.vcf.gz ${CHR} | sort -k 3 > input.vcf
            gsutil -m cp ~{remote_dir}/~{sample_id}/tsvs/clique-representatives-${CHR}.tsv ./representatives.tsv
            bcftools view --header-only all_calls.vcf.gz > ${CHR}.vcf
            join -t $'\t' -1 3 -2 1 input.vcf representatives.tsv | cut -f 1-10 >> ${CHR}.vcf
            bcftools sort ${CHR}.vcf --output-type z > ${CHR}.vcf.gz
            tabix ${CHR}.vcf.gz
            rm -f ${CHR}.vcf
            echo "${CHR}.vcf.gz" >> list.txt
        done
        bcftools concat --file-list list.txt --output-type z > ~{sample_id}.svmerger.vcf.gz
        tabix ~{sample_id}.svmerger.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".svmerger.vcf.gz"
        File output_vcf_gz_tbi = work_dir + "/" + sample_id + ".svmerger.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "8GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}
    