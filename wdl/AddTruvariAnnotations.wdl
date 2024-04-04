version 1.0


# Marks every call as true or false according to $truvari bench$ against a
# truth VCF.
#
workflow AddTruvariAnnotations {
    input {
        String sample_id
        File input_vcf_gz
        File input_vcf_gz_tbi
        File truth_vcf_gz
        File truth_tbi
        String truvari_bench_flags = "--sizemin 0 --sizefilt 0 --sizemax 1000000 --pctsize 0.9 --pctseq 0.9"
    }
    parameter_meta {
    }

    call AddTruvariAnnotationsImpl {
        input:
            sample_id = sample_id,
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            truth_vcf_gz = truth_vcf_gz,
            truth_tbi = truth_tbi,
            truvari_bench_flags = truvari_bench_flags
    }
    
    output {
        File annotated_vcf_gz = AddTruvariAnnotationsImpl.annotated_vcf_gz
        File annotated_tbi = AddTruvariAnnotationsImpl.annotated_tbi
    }
}


task AddTruvariAnnotationsImpl {
    input {
        String sample_id
        File input_vcf_gz
        File input_vcf_gz_tbi
        File truth_vcf_gz
        File truth_tbi
        String truvari_bench_flags
    }
    parameter_meta {
    }
    
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    Int disk_size_gb = 50*( ceil(size(input_vcf_gz,"GB")) + ceil(size(truth_vcf_gz,"GB")) )

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} truvari bench ~{truvari_bench_flags} --base ~{truth_vcf_gz} --comp ~{input_vcf_gz} --output ./truvari/
        
        # Marking TPs
        bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ./truvari/tp-comp.vcf.gz > tmp.tsv
        N_ROWS=$(wc -l < tmp.tsv)
        rm -f bits.txt
        for i in $(seq 1 ${N_ROWS}); do
            echo "1" >> bits.txt
        done
        paste tmp.tsv bits.txt > tps.tsv
        
        # Marking FPs
        bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ./truvari/fp.vcf.gz > tmp.tsv
        N_ROWS=$(wc -l < tmp.tsv)
        rm -f bits.txt
        for i in $(seq 1 ${N_ROWS}); do
            echo "0" >> bits.txt
        done
        paste tmp.tsv bits.txt > fps.tsv
        
        # Annotating
        echo '##INFO=<ID=TRUVARI_TRUE,Number=1,Type=Integer,Description="True according to truvari bench">' > header.txt
        cat tps.tsv fps.tsv | sort -k 1V -k 2n | bgzip -c > annotations.tsv.gz
        bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,INFO/TRUVARI_TRUE ~{input_vcf_gz} --output-type z > ~{sample_id}_annotated.vcf.gz
        tabix -f ~{sample_id}_annotated.vcf.gz
    >>>

    output {
        File annotated_vcf_gz = work_dir + "/" + sample_id + "_annotated.vcf.gz"
        File annotated_tbi = work_dir + "/" + sample_id + "_annotated.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}