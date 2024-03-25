version 1.0


# Runs a truvari inter-sample merge in parallel on every chromosome.
#
workflow TruvariIntersample {
    input {
        String source_dir
        Array[String] chromosomes
        File reference_fai
        File density_counter_py
        Int max_records_per_chunk = 10000
        String destination_dir
        Int monitor_every_seconds = 300
    }
    parameter_meta {
        source_dir: "Contains per-chromosome files built by workflow $FilterAndSplit$."
        destination_dir: "The merged VCFs (one per chromosome) are stored in this remote directory."
        max_records_per_chunk: "Discards chunks that contain more than this many records. Setting it to 10k keeps 99.9% of all chunks in AoU Phase 1 (1027 samples) on CHM13."
        monitor_every_seconds: "Print progress every X seconds"
    }

    scatter (chr in chromosomes) {
        call TruvariIntersampleImpl {
            input:
                source_dir = source_dir,
                chromosome = chr,
                reference_fai = reference_fai,
                density_counter_py = density_counter_py,
                max_records_per_chunk = max_records_per_chunk,
                destination_dir = destination_dir,
                monitor_every_seconds = monitor_every_seconds
        }
    }
    
    output {
    }
}


task TruvariIntersampleImpl {
    input {
        String source_dir
        String chromosome
        File reference_fai
        File density_counter_py
        Int max_records_per_chunk
        String destination_dir
        Int monitor_every_seconds
    }
    parameter_meta {
    }
    
    Int ram_size_gb = 128  # Arbitrary
    String work_dir = "/cromwell_root/truvari_intrasample"
    
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
            TEST=$(gsutil -m cp ~{source_dir}'/*_'~{chromosome}_split.vcf.gz'*' . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        find . -maxdepth 1 -name '*.vcf.gz' > list.txt
        
        # BCFTOOLS
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --force-samples --merge none --file-list list.txt --output-type z > ~{chromosome}.merged.vcf.gz
        tabix -f ~{chromosome}.merged.vcf.gz
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --do-not-normalize --multiallelics -any --output-type z ~{chromosome}.merged.vcf.gz > ~{chromosome}.normed.vcf.gz
        tabix -f ~{chromosome}.normed.vcf.gz

        # Discarding chunks with too many records
        python ~{density_counter_py} ~{chromosome}.normed.vcf.gz > chunks.bed
        awk '$4 >= ~{max_records_per_chunk}' chunks.bed > excluded.bed
        bedtools complement -i excluded.bed -g ~{reference_fai} > included.bed
        bedtools sort -faidx ~{reference_fai} -i included.bed > included.sorted.bed
        
        # TRUVARI
        truvari collapse --input ~{chromosome}.normed.vcf.gz --collapsed-output removed.vcf.gz --sizemin 0 --sizemax 1000000 --keep common --bed included.sorted.bed --gt all \
            | bcftools sort --max-mem $(( ~{ram_size_gb} - 4 ))G --output-type z > ~{chromosome}.collapsed.vcf.gz &
        pid=$!
        while ps -p "${pid}" > /dev/null ; do 
            sleep ~{monitor_every_seconds}
            date
            ls
            bcftools query -f "%CHROM\n" removed.vcf.gz | uniq -c 
        done
        wait
        tabix -f ~{chromosome}.collapsed.vcf.gz
        
        # Uploading
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{chromosome}.collapsed.vcf.'gz*' ~{destination_dir} && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp chunks.bed ~{destination_dir}/~{chromosome}.chunks.bed && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading files. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 8
        memory: ram_size_gb + "GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}