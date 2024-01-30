version 1.0

#
workflow Jasmine {
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
    
    call JasmineImpl {
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
        File output_vcf_gz = JasmineImpl.output_vcf_gz
        File output_vcf_gz_tbi = JasmineImpl.output_vcf_gz_tbi
    }
}


task JasmineImpl {
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
        
        # Removing BND records, since they seem to be a runtime bottleneck.
        touch ~{pbsv_vcf_gz_tbi} ~{sniffles_vcf_gz_tbi} ~{pav_vcf_gz_tbi}
        bcftools view --include 'SVTYPE!="BND"' --output-type z ~{pbsv_vcf_gz} > pbsv.vcf.gz
        tabix pbsv.vcf.gz
        bcftools view --include 'SVTYPE!="BND"' --output-type z ~{sniffles_vcf_gz} > sniffles.vcf.gz
        tabix sniffles.vcf.gz
        bcftools view --include 'SVTYPE!="BND"' --output-type z ~{pav_vcf_gz} > pav.vcf.gz
        tabix pav.vcf.gz
        echo "pbsv.vcf.gz" >> list.txt
        echo "sniffles.vcf.gz" >> list.txt
        echo "pav.vcf.gz" >> list.txt
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list.txt --output-type v > concat.vcf
        
        # Running jasmine with default parameters
        echo "concat.vcf" > list.txt
        chmod +x /opt/conda/envs/jasmine/bin/jasmine
        source activate jasmine
        ${TIME_COMMAND} jasmine --allow_intrasample --output_genotypes threads=${N_THREADS} file_list=list.txt out_file=~{sample_id}.jasmine.vcf
        conda deactivate
        bcftools view --header-only ~{sample_id}.jasmine.vcf
        bcftools sort --output-type z ~{sample_id}.jasmine.vcf > ~{sample_id}.jasmine.vcf.gz
        tabix ~{sample_id}.jasmine.vcf.gz
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".jasmine.vcf.gz"
        File output_vcf_gz_tbi = work_dir + "/" + sample_id + ".jasmine.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: "16GB"  # Arbitrary
        disks: "local-disk 20 HDD"  # Arbitrary
        preemptible: 0
    }
}
