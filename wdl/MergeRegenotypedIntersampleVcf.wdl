version 1.0


# Assume that we have built an inter-sample-merged VCF, that we have
# subsequently kept just one column of it (setting it to all 0/1), and that we
# have re-genotyped this single-column VCF using the BAM of every original
# sample. The program combines all such re-genotyped single-sample VCFs into
# an inter-sample-merged VCF with updated genotypes.
#
# Remark: every file in $regenotyped_vcfs_list$ is assumed to have exactly the
# same set of calls in the same order. This is usually the case if the files
# are the result of re-genotyping the same inter-sample VCF with different BAMs.
# If this is not the case the program gives a wrong output, and it should be
# replaced with $bcftools merge --merge none$.
#
workflow MergeRegenotypedIntersampleVcf {
    input {
        File regenotyped_vcfs_list
    }
    parameter_meta {
    }
    
    call MergeImpl {
        input:
            regenotyped_vcfs_list = regenotyped_vcfs_list
    }
    
    output {
        File output_vcf_gz = MergeImpl.output_vcf_gz
        File output_tbi = MergeImpl.output_tbi
    }
}


task MergeImpl {
    input {
        File regenotyped_vcfs_list
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$((2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Initializing the inter-sample VCF with the first file
        ADDRESS=$(head -n 1 ~{regenotyped_vcfs_list})
        gsutil -m cp ${ADDRESS} first.vcf.gz
        bcftools view --header-only first.vcf.gz > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > header.txt
        FIELDS=$(tail -n 1 tmp.txt)
        bcftools view --no-header first.vcf.gz > body.txt
        rm -f first.vcf.gz
        
        # Appending all remaining files
        N_FILES=$(wc -l < ~{regenotyped_vcfs_list})
        tail -n $(( ${N_FILES} - 1 )) ~{regenotyped_vcfs_list} > list.txt
        while read ADDRESS; do
            # Adding the new sample to the set of columns
            gsutil -m cp ${ADDRESS} ./tmp.vcf.gz
            bcftools view --header-only tmp.vcf.gz > tmp.txt
            N_ROWS=$(wc -l < tmp.txt)
            SAMPLE_ID=$(tail -n 1 tmp.txt | cut -f 10)
            FIELDS=$(echo -e "${FIELDS}\t${SAMPLE_ID}")
            echo "Current columns: ${FIELDS}"
            rm -f tmp.txt
            # Adding the new column to the body
            bcftools view --no-header tmp.vcf.gz | cut -f 10 > column.txt
            paste body.txt column.txt > body_new.txt
            rm -f body.txt; mv body_new.txt body.txt
            rm -f tmp.vcf.gz
        done < list.txt
        
        # Merging
        echo -e "${FIELDS}" > fields.txt
        cat header.txt fields.txt body.txt > merged.vcf
        rm -f *.txt
        ${TIME_COMMAND} bgzip -@ ${N_THREADS} merged.vcf
        tabix -f merged.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 8
        memory: "32GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}