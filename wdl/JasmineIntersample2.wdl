version 1.0

#
workflow JasmineIntersample2 {
    input {
        Array[File] input_vcf_gz
        Array[File] input_tbi
    }
    parameter_meta {
    }
    
    Array[String] chromosomes = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
    scatter(chr in chromosomes) {
        call JasmineIntersample2Impl {
            input:
                input_vcf_gz = input_vcf_gz,
                input_tbi = input_tbi,
                chromosome = chr
        }
    }
    call ConcatChromosomeVCFs {
        input:
            chromosome_vcf_gz = JasmineIntersample2Impl.output_vcf_gz,
            chromosome_tbi = JasmineIntersample2Impl.output_tbi
    }
    
    output {
        File output_vcf_gz = ConcatChromosomeVCFs.output_vcf_gz
        File output_tbi = ConcatChromosomeVCFs.output_tbi
    }
}


task JasmineIntersample2Impl {
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
        
        # Removing BND records, since they seem to be a runtime bottleneck.
        INPUT_FILES=~{sep=',' input_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            ID=$(basename ${INPUT_FILE} .vcf.gz)
            gunzip -c ${INPUT_FILE} > ${ID}.vcf
            tail -n 10 ${ID}.vcf
            bgzip ${ID}.vcf
            tabix ${ID}.vcf.gz
            bcftools view --include 'SVTYPE!="BND"' ${ID}.vcf.gz ~{chromosome} > ${ID}.filtered.vcf
            tail -n 10 ${ID}.filtered.vcf
            echo ${ID}.filtered.vcf >> list.txt
            rm -f ${ID}.vcf*
        done
        
        # Running jasmine with default parameters
        chmod +x /opt/conda/envs/jasmine/bin/jasmine
        source activate jasmine
        ${TIME_COMMAND} jasmine --output_genotypes threads=${N_THREADS} file_list=list.txt out_file=~{chromosome}.vcf
        conda deactivate
        bcftools sort --output-type z ~{chromosome}.vcf > ~{chromosome}.vcf.gz
        tabix ~{chromosome}.vcf.gz
    >>>
    output {
        File output_vcf_gz = work_dir + "/" + chromosome + ".vcf.gz"
        File output_tbi = work_dir + "/" + chromosome + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "64GB"  # Arbitrary
        disks: "local-disk 500 HDD"  # Arbitrary
        preemptible: 0
    }
}


# Creates a single merged VCF from the concatenation of all single-chromosome
# VCFs produced by jasmine.
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
