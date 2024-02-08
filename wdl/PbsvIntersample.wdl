version 1.0


#
workflow PbsvIntersample {
    input {
        Array[File] input_bam
        Array[File] input_bai
        Array[String] input_ids
        File tandems_bed
        File reference_fa
    }

    scatter (i in range(length(input_bam))) {
        call SingleSampleCalling {
            input:
                input_bam = input_bam[i],
                input_bai = input_bai[i],
                sample_id = input_ids[i],
                tandems_bed = tandems_bed
        }
    }
    call JointCalling {
        input:
            input_svsig = SingleSampleCalling.svsig,
            reference_fa = reference_fa
    }

    output {
         File output_vcf_gz = JointCalling.output_vcf_gz
         File output_tbi = JointCalling.output_tbi
    }
}


task SingleSampleCalling {
    input {
        File input_bam
        File input_bai
        String sample_id
        File tandems_bed
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

        ${TIME_COMMAND} pbsv discover --ccs --tandem-repeats ~{tandems_bed} ~{input_bam} ~{sample_id}.svsig.gz
        ls -laht
        tree
    >>>

    output {
        File svsig = work_dir + "/" + sample_id + ".svsig.gz"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "32GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0
    }
}


task JointCalling {
    input {
        Array[File] input_svsig
        File reference_fa
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
        
        INPUT_FILES=~{sep=',' input_svsig}
        INPUT_FILES=$(echo ${INPUT_FILES} | tr ',' ' ')
        ${TIME_COMMAND} pbsv call --num-threads ${N_THREADS} --ccs ~{reference_fa} ${INPUT_FILES} joint.vcf
        bgzip joint.vcf
        tabix joint.vcf.gz
        ls -laht
        tree
    >>>

    output {
        File output_vcf_gz = work_dir + "/joint.vcf.gz"
        File output_tbi = work_dir + "/joint.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 32
        memory: "64GB"  # Arbitrary
        disks: "local-disk 50 HDD"  # Arbitrary
        preemptible: 0
    }
}