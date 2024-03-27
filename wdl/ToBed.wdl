version 1.0


#
workflow ToBed {
    input {
        File vcf_gz_file
        File tbi_file
        File bed_file
    }
    parameter_meta {
    }

    call ToBedImpl {
        input:
            vcf_gz_file = vcf_gz_file,
            tbi_file = tbi_file,
            bed_file = bed_file
    }
    
    output {
        File vcf_gz = ToBedImpl.vcf_gz
        File vcf_gz_tbi = ToBedImpl.vcf_gz_tbi
    }
}


task ToBedImpl {
    input {
        File vcf_gz_file
        File tbi_file
        File bed_file
    }
    parameter_meta {
    }
    
    String docker_dir = "/truvari_intrasample"
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
        
        bcftools view --header-only ~{vcf_gz_file} > out.vcf
        bedtools intersect -u -a ~{vcf_gz_file} -b ~{bed_file} >> out.vcf
        bgzip out.vcf
        tabix -f out.vcf.gz
    >>>

    output {
        File vcf_gz = work_dir + "/out.vcf.gz"
        File vcf_gz_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 1
        memory: "16GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}