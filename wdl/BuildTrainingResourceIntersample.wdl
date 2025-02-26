version 1.0


#
workflow BuildTrainingResourceIntersample {
    input {
        Array[File] intrasample_vcf_gz
        Array[File] intrasample_tbi
    }
    parameter_meta {
    }
    
    call BuildTrainingResourceIntersampleImpl {
        input:
            intrasample_vcf_gz = intrasample_vcf_gz,
            intrasample_tbi = intrasample_tbi
    }
    
    output {
        File merged_vcf = BuildTrainingResourceIntersampleImpl.merged_vcf
        File merged_tbi = BuildTrainingResourceIntersampleImpl.merged_tbi
    }
}


#
task BuildTrainingResourceIntersampleImpl {
    input {
        Array[File] intrasample_vcf_gz
        Array[File] intrasample_tbi
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        INPUT_FILES=~{sep=',' intrasample_vcf_gz}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        rm -f list.txt
        for INPUT_FILE in ${INPUT_FILES}; do
            echo ${INPUT_FILE} >> list.txt
        done
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list.txt --output-type z > merged.vcf.gz
        tabix -f merged.vcf.gz
    >>>

    output {
        File merged_vcf = work_dir + "/merged.vcf.gz"
        File merged_tbi = work_dir + "/merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "32GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
