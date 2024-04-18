version 1.0


workflow Hapdiff {
    input {
        File ref_fa
        File tandems_bed = "null"
        File pat_fa
        File mat_fa
        Int min_sv_length = 50
        Int fragment_size_mb = 0
        Int n_cpu = 32
        Int mem_gb = 128
    }
    parameter_meta {
        fragment_size_mb: "0: do not fragment"
        tandems_bed: "null: do not use"
    }

    call HapdiffImpl {
        input:
            ref_fa=ref_fa,
            tandems_bed=tandems_bed,
            pat_fa=pat_fa,
            mat_fa=mat_fa,
            min_sv_length=min_sv_length,
            fragment_size_mb=fragment_size_mb,
            n_cpu=n_cpu,
            mem_gb=mem_gb
    }

    output {
        File phased_vcf = HapdiffImpl.phased_vcf
        File phased_tbi = HapdiffImpl.phased_tbi
        File confident_bed = HapdiffImpl.confident_bed
    }
}


task HapdiffImpl {
    input {
        File ref_fa
        File tandems_bed
        File pat_fa
        File mat_fa
        Int min_sv_length
        Int fragment_size_mb
        Int n_cpu
        Int mem_gb
    }
    
    String docker_dir = "/hapdiff"
    String work_dir = "/cromwell_root/truvari_intrasample"
    String output_dir = "/cromwell_root/truvari_intrasample/output"
    Int disk_size_gb = 20*( ceil(size(ref_fa,"GB")) + ceil(size(pat_fa,"GB")) + ceil(size(mat_fa,"GB")) )

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        if [ ~{tandems_bed} = "null" ]; then
             TANDEMS_FLAG=""
        else
             TANDEMS_FLAG="--tandem-repeats ~{tandems_bed}"
        fi
        if [ ~{fragment_size_mb} = "0" ]; then
             FRAGMENT_FLAG=""
        else
             FRAGMENT_FLAG="--fragment ~{fragment_size_mb}"
        fi
        ${TIME_COMMAND} python3 ~{docker_dir}/hapdiff.py \
            --reference ~{ref_fa} \
            --pat ~{pat_fa} \
            --mat ~{mat_fa} \
            --out-dir ~{output_dir} \
            ${TANDEMS_FLAG} \
            --sv-size ~{min_sv_length} \
            ${FRAGMENT_FLAG} \
            --threads ${N_THREADS}
        tabix -f ${OUTPUT_DIR}/hapdiff_phased.vcf.gz
    >>>
    
    output {
        File phased_vcf = output_dir + "/hapdiff_phased.vcf.gz"
        File phased_tbi = output_dir + "/hapdiff_phased.vcf.gz.tbi"
        File confident_bed = output_dir + "/confident_regions.bed"
    }
    runtime {
        docker: "mkolmogo/hapdiff:0.9"
        cpu: n_cpu
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
    }
}
