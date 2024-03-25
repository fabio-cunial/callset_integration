version 1.0


# Applies a filter to a VCF and sequentially splits it into a given set of 
# chromosomes.
#
task FilterAndSplit {
    input {
        String sample_id
        File sample_vcf_gz
        File sample_vcf_gz_tbi
        String filter_string
        Array[String] chromosomes
        String destination_dir
    }
    parameter_meta {
        sample_vcf_gz: "Assumed to be already sorted"
        filter_string: "String to be used in $bcftools filter --include$. Empty = no filtering."
        destination_dir: "The filtered and split files are stored in this remote directory."
    }
    
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    Int disk_size_gb = 10*(ceil(size(sample_vcf_gz,"GB")))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
                
        CHROMOSOMES=~{sep='-' chromosomes}
        CHROMOSOMES=$(echo ${CHROMOSOMES} | tr '-' ' ')
        if [ -z ~{filter_string} ]; then
            INCLUDE_STR=""
        else
            INCLUDE_STR="--include ~{filter_string}"
        fi
        for CHROMOSOME in ${CHROMOSOMES}; do
            ${TIME_COMMAND} bcftools filter ${INCLUDE_STR} --regions ${CHROMOSOME} --output-type z ~{sample_vcf_gz} > ~{sample_id}_${CHROMOSOME}_split.vcf.gz
            tabix -f ~{sample_id}_${CHROMOSOME}_split.vcf.gz
        done
        gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp '*_split.vcf.gz*' ~{destination_dir}
    >>>

    output {
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}