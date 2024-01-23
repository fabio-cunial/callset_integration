version 1.0


#
workflow Svimmer {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
    }
    parameter_meta {
    }
    call SvimmerImpl {
        input:
            sample_id = sample_id,
            pbsv_vcf_gz = pbsv_vcf_gz,
            pbsv_vcf_gz_tbi = pbsv_vcf_gz_tbi,
            sniffles_vcf_gz = sniffles_vcf_gz,
            sniffles_vcf_gz_tbi = sniffles_vcf_gz_tbi,
            pav_vcf_gz = pav_vcf_gz,
            pav_vcf_gz_tbi = pav_vcf_gz_tbi
    }
    output {
        File output_vcf_gz = SvimmerImpl.output_vcf_gz
        File output_vcf_gz_tbi = SvimmerImpl.output_vcf_gz_tbi
    }
}


task SvimmerImpl {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
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
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        TIME_COMMAND="/usr/bin/time --verbose"
        
        mv ~{pbsv_vcf_gz} . ; mv ~{pbsv_vcf_gz_tbi} .
        mv ~{sniffles_vcf_gz} . ; mv ~{sniffles_vcf_gz_tbi} .
        mv ~{pav_vcf_gz} . ; mv ~{pav_vcf_gz_tbi} .
        ls *.vcf.gz > list.txt
        cat list.txt
        touch *.tbi
        ${TIME_COMMAND} svimmer --threads ${N_THREADS} --ids --output ~{sample_id}.svimmer.vcf list.txt chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
        ls -laht
        bcftools sort ~{sample_id}.svimmer.vcf --output-type z > ~{sample_id}.svimmer.vcf.gz
        tabix ~{sample_id}.svimmer.vcf.gz        
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".svimmer.vcf.gz"
        File output_vcf_gz_tbi = work_dir + "/" + sample_id + ".svimmer.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "64GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0
    }
}
