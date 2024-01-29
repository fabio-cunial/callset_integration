version 1.0


#
workflow Svimmer {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
    }
    parameter_meta {
    }
    Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    scatter(chr in chromosomes) {
        call SvimmerImpl {
            input:
                input_vcf_gz = input_vcf_gz,
                input_tbi = input_tbi,
                chromosome = chr
        }
    }
    call ConcatChromosomeVCFs {
        input:
            chromosome_vcf_gz = SvimmerImpl.output_vcf_gz,
            chromosome_tbi = SvimmerImpl.output_tbi
    }
    
    output {
        File output_vcf_gz = ConcatChromosomeVCFs.output_vcf_gz
        File output_vcf_gz_tbi = ConcatChromosomeVCFs.output_tbi
    }
}


task SvimmerImpl {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
        String chromosome
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
        
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} python3 ~{docker_dir}/svimmer/svimmer --threads ${N_THREADS} --ids --output ~{chromosome}.vcf list.txt ~{chromosome}
        bgzip ~{chromosome}.vcf
        tabix ~{chromosome}.vcf.gz        
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + chromosome + ".vcf.gz"
        File output_tbi = work_dir + "/" + chromosome + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: "32GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0
    }
}


# Creates a single merged VCF from the concatenation of all single-chromosome
# VCFs produced by svimmer.
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
        i="0"
        for INPUT_FILE in ${INPUT_FILES}; do
            bcftools view --header-only ${INPUT_FILE} > header1.txt
            N_ROWS=$(wc -l < header1.txt)
            head -n $(( ${N_ROWS} - 1 )) header1.txt > fixed${i}.vcf
            echo "##INFO=<ID=ID,Number=1,Type=String,Description="id">" >> fixed${i}.vcf
            echo "##INFO=<ID=TIG_REGION,Number=1,Type=String,Description="tig">" >> fixed${i}.vcf
            echo "##INFO=<ID=QUERY_STRAND,Number=1,Type=String,Description="strand">" >> fixed${i}.vcf
            echo "##INFO=<ID=HOM_REF,Number=1,Type=String,Description="homref">" >> fixed${i}.vcf
            echo "##INFO=<ID=HOM_TIG,Number=1,Type=String,Description="homtig">" >> fixed${i}.vcf
            echo "##INFO=<ID=PRECISE,Number=1,Type=String,Description="precise">" >> fixed${i}.vcf
            echo "##INFO=<ID=SUPPORT,Number=1,Type=Integer,Description="support">" >> fixed${i}.vcf
            echo "##INFO=<ID=COVERAGE,Number=1,Type=String,Description="coverage">" >> fixed${i}.vcf
            echo "##INFO=<ID=STRAND,Number=1,Type=String,Description="strand">" >> fixed${i}.vcf
            echo "##INFO=<ID=AF,Number=1,Type=Float,Description="af">" >> fixed${i}.vcf
            echo "##INFO=<ID=STDEV_LEN,Number=1,Type=Float,Description="std">" >> fixed${i}.vcf
            echo "##INFO=<ID=STDEV_POS,Number=1,Type=Float,Description="std">" >> fixed${i}.vcf
            echo "##INFO=<ID=SUPPORT_LONG,Number=1,Type=Integer,Description="support">" >> fixed${i}.vcf
            echo "##INFO=<ID=CHR2,Number=1,Type=String,Description="chr2">" >> fixed${i}.vcf
            echo "##INFO=<ID=INNER_REF,Number=1,Type=String,Description="inner">" >> fixed${i}.vcf
            echo "##INFO=<ID=INNER_TIG,Number=1,Type=String,Description="inner">" >> fixed${i}.vcf
            tail -n 1 header1.txt >> fixed${i}.vcf
            bcftools view --no-header ${INPUT_FILE} >> fixed${i}.vcf
            bgzip fixed${i}.vcf
            tabix fixed${i}.vcf.gz
            echo fixed${i}.vcf.gz >> list.txt
            i=$(( ${i} + 1 ))
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
