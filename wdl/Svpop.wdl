version 1.0


#
workflow Svpop {
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


task SvpopImpl {
    input {
        String sample_id
        File pbsv_vcf_gz
        File pbsv_vcf_gz_tbi
        File sniffles_vcf_gz
        File sniffles_vcf_gz_tbi
        File pav_vcf_gz
        File pav_vcf_gz_tbi
        File reference_fa
        File reference_fai
        String ucsc_ref_name
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
        
        touch ~{pbsv_vcf_gz_tbi} ~{sniffles_vcf_gz_tbi} ~{pav_vcf_gz_tbi}        
        
        
        
        # Parameters are set to approximate truvari's defaults.
        rm -rf config/; mkdir config/
        touch config/samples.tsv
        echo -e "NAME\tSAMPLE\tTYPE\tDATA\tVERSION\tPARAMS\tCOMMENT" >> config/samples.tsv
        echo -e "pbsv\tDEFAULT\tpbsv\t~{pbsv_vcf_gz}\t1\t\t"  >> config/samples.tsv
        echo -e "sniffles\tDEFAULT\tsniffles2\t~{sniffles_vcf_gz}\t1\t\t"  >> config/samples.tsv
        echo -e "pav\tDEFAULT\tvcf\t~{pav_vcf_gz}\t1\t\t"  >> config/samples.tsv
        cat config/samples.tsv
        touch config/config.json
        echo "{" >> config/config.json
        echo "\"reference\": \"~{reference_fa}\"," >> config/config.json
        echo "\"reference_fai\": \"~{reference_fai}\"," >> config/config.json
        echo "\"ucsc_ref_name\": \"~{ucsc_ref_name}\"," >> config/config.json
        # samplelist section
        echo "\"samplelist\": {" >> config/config.json
        echo "\"mySamples\": [" >> config/config.json
        i=0
        while read VCF_FILE; do
           if [ $i -eq 0 ]; then
               echo -n "\"${VCF_FILE%.vcf.gz}\"" >> config/config.json
               i=1
           else
               echo -en ",\n\"${VCF_FILE%.vcf.gz}\"" >> config/config.json
           fi
        done < list.txt
        echo -e "\n]" >> config/config.json
        echo "}," >> config/config.json
        # sampleset section
        echo "\"sampleset\": {" >> config/config.json
        echo "\"myMerge\": {" >> config/config.json
        echo "\"sourcetype\": \"caller\"," >> config/config.json
        echo "\"sourcename\": \"sniffles2\"," >> config/config.json
        echo "\"merge\": \"nr::szro(szro=0.7,dist=500,match(score=0.7,limit=4000,ksize=9))\"," >> config/config.json
        echo "\"name\": \"myMerge\"," >> config/config.json
        echo "\"description\": \"myMerge\"" >> config/config.json
        echo "}," >> config/config.json
        echo "}" >> config/config.json
        # end of config file
        echo "}" >> config/config.json
        cat config/config.json
        source activate svpop
        # INS return error at line: <https://github.com/EichlerLab/svpop/blob/1d62b72187172a7898443c5b1a27f3838621a199/svpoplib/svmerge.py#L822>
        ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_ins.bed.gz
        # DEL returns error at line: <https://github.com/EichlerLab/svpop/blob/1d62b72187172a7898443c5b1a27f3838621a199/svpoplib/svmerge.py#L822>
        ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_del.bed.gz
        ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_inv.bed.gz
        ${TIME_COMMAND} snakemake -s ~{docker_dir}/svpop/Snakefile --cores ${N_THREADS} results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_dup.bed.gz
        conda deactivate
        zcat results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_ins.bed.gz \
            results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_del.bed.gz \
            results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_inv.bed.gz \
            results/variant/sampleset/myMerge/mySamples/all/all/bed/sv_dup.bed.gz | gzip > svpop.bed.gz
        uploadVCF svpop.bed.gz " "
        
        
        
        ls -laht
        bgzip ~{sample_id}.svpop.vcf
        tabix ~{sample_id}.svpop.vcf.gz        
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + sample_id + ".svpop.vcf.gz"
        File output_vcf_gz_tbi = work_dir + "/" + sample_id + ".svpop.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: "16GB"  # Arbitrary
        disks: "local-disk 20 HDD"  # Arbitrary
        preemptible: 0
    }
}
