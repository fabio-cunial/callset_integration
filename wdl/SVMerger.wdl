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
    #Array[String] chromosomes = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"]
    Array[String] chromosomes = ["chr1"]
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
    
    output {
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
            grep ^chr${CHR} pbsv.tsv >> chr${CHR}.tsv
            grep ^chr${CHR} sniffles.tsv >> chr${CHR}.tsv
            grep ^chr${CHR} pav.tsv >> chr${CHR}.tsv
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
            
            if [ ${SVTYPE} = "del" ]; then
                AWK_COMMAND='$7 == "DEL"'
            else
                AWK_COMMAND='$7 == "INS"'
            fi
            awk "${AWK_COMMAND}" ${CHROMOSOME_ID}.tsv > ${CHROMOSOME_ID}.${SVTYPE}.tsv
            ${TIME_COMMAND} python2 ~{docker_dir}/sv-merger/main.py MERGE ${CHROMOSOME_ID}.del.tsv ${CHROMOSOME_ID}.trf.sorted.gor ${SVTYPE}
            sort -k 4 ${CHROMOSOME_ID}.${SVTYPE}.tsv > chr.tsv
            ls -laht; tree
            # Outside TRs
            sort -k 1 ${CHROMOSOME_ID}.tsv.outtrr.merged.csv > outtrr.tsv
            join -1 4 -2 1 chr.tsv outtrr.tsv | tr ' ' '\t' | sort --version-sort --key 9 > cliques.tsv
            ls -laht; tree
            java SVMergerGetRepresentative cliques.tsv > clique-representatives-${SVTYPE}-outttr.tsv
            # Inside TRs
            sort -k 1 ${CHROMOSOME_ID}.tsv.intrr.merged.csv > intrr.tsv
            join -1 4 -2 1 chr.tsv intrr.tsv | tr ' ' '\t' | sort --version-sort --key 9 > cliques.tsv
            ls -laht; tree
            java SVMergerGetRepresentative cliques.tsv > clique-representatives-${SVTYPE}-inttr.tsv
            cat clique-representatives-${SVTYPE}-outttr.tsv clique-representatives-${SVTYPE}-inttr.tsv | sort > clique-representatives-${SVTYPE}.tsv
        }
        
        cp ~{docker_dir}/*.class .
        gsutil -m cp ~{remote_dir}/~{sample_id}/tsvs/~{chromosome_id}.tsv .
        gsutil -m cp ~{remote_trf_dir}/~{chromosome_id}.trf.sorted.gor  .        
        svMerger ~{chromosome_id} del
        ls -laht; tree
        svMerger ~{chromosome_id} ins
        ls -laht; tree
        cat clique-representatives-del.tsv clique-representatives-ins.tsv | sort > clique-representatives-~{chromosome_id}.tsv
    >>>
    
    output {
        File outtrr = work_dir + "/clique-representatives-" + chromosome_id + ".tsv"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "8GB"
        disks: "local-disk 50 HDD"
        preemptible: 0
    }
}
