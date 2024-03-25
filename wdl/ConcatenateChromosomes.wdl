version 1.0


# Simply concatenates per-chromosome intersample-merged VCFs into a single VCF.
#
workflow ConcatenateChromosomes {
    input {
        String source_dir
    }
    parameter_meta {
        source_dir: "The output dir of workflow $TruvariIntersample$, containing sorted per-chromosome VCFs."
    }

    call ConcatenateChromosomesImpl {
        input:
            source_dir = source_dir
    }
    
    output {
        File vcf_gz = ConcatenateChromosomesImpl.vcf_gz
        File vcf_gz_tbi = ConcatenateChromosomesImpl.vcf_gz_tbi
    }
}


task ConcatenateChromosomesImpl {
    input {
        String source_dir
    }
    parameter_meta {
    }
    
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    Int disk_size_gb = 256  # Arbitrary

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        GSUTIL_UPLOAD_THRESHOLD="-o GSUtil:parallel_composite_upload_threshold=150M"
        GSUTIL_DELAY_S="600"
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Downloading
        while : ; do
            TEST=$(gsutil -m cp ~{source_dir}'/*.vcf.gz' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        find . -depth 1 -type f -name "*.vcf.gz" > list.txt
        
        ${TIME_COMMAND} bcftools concat --threads ${N_THREADS} --allow-overlaps --file-list list.txt --output-type z > concat.vcf.gz
        tabix -f concat.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/concat.vcf.gz"
        File vcf_gz_tbi = work_dir + "/concat.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 16
        memory: "128GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}