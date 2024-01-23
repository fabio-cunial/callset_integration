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
        ls -laht ~{docker_dir}/svimmer/
        cat ~{docker_dir}/svimmer/svimmer
        python3 ~{docker_dir}/svimmer/svimmer -h || echo 1
        ${TIME_COMMAND} python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer.vcf list.txt chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer1.vcf list.txt chr1
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer2.vcf list.txt chr2
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer3.vcf list.txt chr3
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer4.vcf list.txt chr4
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer5.vcf list.txt chr5
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer6.vcf list.txt chr6
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer7.vcf list.txt chr7
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer8.vcf list.txt chr8
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer9.vcf list.txt chr9
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer10.vcf list.txt chr10
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer11.vcf list.txt chr11
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer12.vcf list.txt chr12
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer13.vcf list.txt chr13
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer14.vcf list.txt chr14
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer15.vcf list.txt chr15
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer16.vcf list.txt chr16
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer17.vcf list.txt chr17
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer18.vcf list.txt chr18
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer19.vcf list.txt chr19
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer20.vcf list.txt chr20
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer21.vcf list.txt chr21
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmer22.vcf list.txt chr22
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmerX.vcf list.txt chrX
        python3 ~{docker_dir}/svimmer/svimmer --ids --output ~{sample_id}.svimmerY.vcf list.txt chrY
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
