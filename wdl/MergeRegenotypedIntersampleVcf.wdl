version 1.0


# Assume that we have built an inter-sample-merged VCF, that we have
# subsequently kept just one column of it (setting it to all 0/1), and that we
# have re-genotyped this single-column VCF using the BAM of every original
# sample. The program combines all such re-genotyped single-sample VCFs into
# an inter-sample-merged VCF with updated genotypes. 
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


#
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
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        # Downloading
        rm -f list.txt
        i="0"
        while read ADDRESS; do
            i=$(( ${i} + 1 ))
            gsutil -m cp ${ADDRESS} ./${i}.vcf.gz
            tabix -f ${i}.vcf.gz
            echo ${i}.vcf.gz >> list.txt
        done < ~{regenotyped_vcfs_list}
        
        # Merging
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --file-list list.txt --output-type z > merged.vcf.gz
        tabix -f merged.vcf.gz
    >>>

    output {
        File output_vcf_gz = work_dir + "/merged.vcf.gz"
        File output_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 16
        memory: "64GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}