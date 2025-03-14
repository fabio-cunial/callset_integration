version 1.0


#
workflow Phase2RocInvestigation {
    input {
        Int min_sv_length
        
        Array[String] sample_id
        Array[File] intrasample_vcf
        Array[File] intrasample_tbi
        Array[File] dipcall_asm10_vcf
        Array[File] dipcall_asm10_tbi
        
        File dipcall_asm10_intersample_vcf
        File dipcall_asm10_intersample_tbi
        
        File truth_bed
        
        File training_resource_12x_vcf
        File training_resource_12x_tbi
        File training_resource_25x_vcf
        File training_resource_25x_tbi
        File training_resource_old_08x_vcf
        File training_resource_old_08x_tbi
        
        File dipcall_old_intersample_vcf
        File dipcall_old_intersample_tbi
        
    }
    parameter_meta {
    }
    
    
    scatter(i in range(length(intrasample_vcf))) {
        call CountCalls as count1 {
            input: 
                min_sv_length = min_sv_length,
                sample_id = sample_id[i],
                single_sample_vcf = intrasample_vcf[i],
                single_sample_tbi = intrasample_tbi[i]
        }
        call BenchSingleSample as bench1 {
            input: 
                min_sv_length = min_sv_length,
                sample_id = sample_id[i],
                single_sample_vcf = intrasample_vcf[i],
                single_sample_tbi = intrasample_tbi[i],
                dipcall_asm10_vcf = dipcall_asm10_vcf[i],
                dipcall_asm10_tbi = dipcall_asm10_tbi[i],
                dipcall_asm10_intersample_vcf = dipcall_asm10_intersample_vcf,
                dipcall_asm10_intersample_tbi = dipcall_asm10_intersample_tbi,
                truth_bed = truth_bed
        }
    }
    call CountCalls as count2 {
        input: 
            min_sv_length = min_sv_length,
            sample_id = "training_resource_12x",
            single_sample_vcf = training_resource_12x_vcf,
            single_sample_tbi = training_resource_12x_tbi
    }
    call CountCalls as count3 {
        input: 
            min_sv_length = min_sv_length,
            sample_id = "training_resource_25x",
            single_sample_vcf = training_resource_25x_vcf,
            single_sample_tbi = training_resource_25x_tbi
    }
    call CountCalls as count4 {
        input: 
            min_sv_length = min_sv_length,
            sample_id = "training_resource_old_08x",
            single_sample_vcf = training_resource_old_08x_vcf,
            single_sample_tbi = training_resource_old_08x_tbi
    }
    call BenchTrainingResources as bench2 {
        input:
            min_sv_length = min_sv_length,
            training_resource_12x_vcf = training_resource_12x_vcf,
            training_resource_12x_tbi = training_resource_12x_tbi,
            training_resource_25x_vcf = training_resource_25x_vcf,
            training_resource_25x_tbi = training_resource_25x_tbi,
            training_resource_old_08x_vcf = training_resource_old_08x_vcf,
            training_resource_old_08x_tbi = training_resource_old_08x_tbi,
            truth_bed = truth_bed
    }
    call CountCalls as count5 {
        input: 
            min_sv_length = min_sv_length,
            sample_id = "dipcall_old_intersample",
            single_sample_vcf = dipcall_old_intersample_vcf,
            single_sample_tbi = dipcall_old_intersample_tbi
    }
    call CountCalls as count6 {
        input: 
            min_sv_length = min_sv_length,
            sample_id = "dipcall_asm10_intersample",
            single_sample_vcf = dipcall_asm10_intersample_vcf,
            single_sample_tbi = dipcall_asm10_intersample_tbi
    }
    call BenchDipcalls as bench3 {
        input:
            dipcall_old_intersample_vcf = dipcall_old_intersample_vcf,
            dipcall_old_intersample_tbi = dipcall_old_intersample_tbi,
            dipcall_asm10_intersample_vcf = dipcall_asm10_intersample_vcf,
            dipcall_asm10_intersample_tbi = dipcall_asm10_intersample_tbi,
            truth_bed = truth_bed
    }
    
    output {
    }
}


#
task CountCalls {
    input {
        Int min_sv_length
        
        String sample_id
        File single_sample_vcf
        File single_sample_tbi
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*ceil(size(single_sample_vcf, "GB")) + 50
    String prefix = "~{sample_id}_~{min_sv_length}"
    String n_ins = "~{prefix}_ins.txt"
    String n_del = "~{prefix}_del.txt"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVTYPE=\"INS\" && (SVLEN<=-~{min_sv_length} || SVLEN>=~{min_sv_length})" --output-type z ~{single_sample_vcf} > tmp.vcf.gz
        tabix -f tmp.vcf.gz
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --no-header tmp.vcf.gz | wc -l > ~{n_ins}
        ${TIME_COMMAND} bcftools filter --threads ${N_THREADS} --include "SVTYPE=\"DEL\" && (SVLEN<=-~{min_sv_length} || SVLEN>=~{min_sv_length})" --output-type z ~{single_sample_vcf} > tmp.vcf.gz
        tabix -f tmp.vcf.gz
        ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --no-header tmp.vcf.gz | wc -l > ~{n_del}
    >>>

    output {
        File n_ins = work_dir + "/" + n_ins
        File n_del = work_dir + "/" + n_del
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: "16GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task BenchSingleSample {
    input {
        Int min_sv_length
        
        String sample_id
        File single_sample_vcf
        File single_sample_tbi
        
        File dipcall_asm10_vcf
        File dipcall_asm10_tbi
        File dipcall_asm10_intersample_vcf
        File dipcall_asm10_intersample_tbi 
        
        File truth_bed
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*ceil(size(single_sample_vcf, "GB") + size(dipcall_asm10_intersample_vcf, "GB")) + 50
    Int ram_gb = 32
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        TRUVARI_BENCH_SETTINGS_1="--sizefilt 0 --sizemin 0 --sizemax 1000000"
        TRUVARI_BENCH_SETTINGS_2="--sizefilt ~{min_sv_length} --sizemin ~{min_sv_length} --sizemax 1000000 --pctsize 0.9 --pctseq 0.9"
        
        mv ~{single_sample_vcf} calls.vcf.gz
        mv ~{single_sample_tbi} calls.vcf.gz.tbi
        
        # Intrasample
        ${TIME_COMMAND} bcftools norm --threads ${N_THREADS} --multiallelics - --output-type z ~{dipcall_asm10_vcf} > truth.vcf.gz
        tabix -f truth.vcf.gz
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_1} --includebed ~{truth_bed} -b truth.vcf.gz -c calls.vcf.gz -o truvari/
        mv truvari/summary.json ./~{sample_id}_1.json
        rm -rf truvari/
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_2} --includebed ~{truth_bed} -b truth.vcf.gz -c calls.vcf.gz -o truvari/
        mv truvari/summary.json ./~{sample_id}_2.json
        rm -rf truvari/
        
        # Intersample
        mv ~{dipcall_asm10_intersample_vcf} truth.vcf.gz
        mv ~{dipcall_asm10_intersample_tbi} truth.vcf.gz.tbi
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_1} --includebed ~{truth_bed} -b truth.vcf.gz -c calls.vcf.gz -o truvari/
        mv truvari/summary.json ./~{sample_id}_intersample_1.json
        rm -rf truvari/
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_2} --includebed ~{truth_bed} -b truth.vcf.gz -c calls.vcf.gz -o truvari/
        mv truvari/summary.json ./~{sample_id}_intersample_2.json
        rm -rf truvari/
    >>>

    output {
        File summary_1 = work_dir + "/" + sample_id + "_1.json"
        File summary_2 = work_dir + "/" + sample_id + "_2.json"
        File summary_intersample_1 = work_dir + "/" + sample_id + "_intersample_1.json"
        File summary_intersample_2 = work_dir + "/" + sample_id + "_intersample_2.json"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task BenchTrainingResources {
    input {
        Int min_sv_length
        
        File training_resource_12x_vcf
        File training_resource_12x_tbi
        File training_resource_25x_vcf
        File training_resource_25x_tbi
        File training_resource_old_08x_vcf
        File training_resource_old_08x_tbi
        
        File truth_bed
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*ceil(size(training_resource_12x_vcf, "GB") + size(training_resource_25x_vcf, "GB") + size(training_resource_old_08x_vcf, "GB")) + 50
    Int ram_gb = 32
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        TRUVARI_BENCH_SETTINGS_1="--sizefilt 0 --sizemin 0 --sizemax 1000000"
        TRUVARI_BENCH_SETTINGS_2="--sizefilt ~{min_sv_length} --sizemin ~{min_sv_length} --sizemax 1000000 --pctsize 0.9 --pctseq 0.9"
        
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_1} --includebed ~{truth_bed} -b ~{training_resource_old_08x_vcf} -c ~{training_resource_12x_vcf} -o truvari/
        mv truvari/summary.json ./12x_vs_8x_1.json
        rm -rf truvari/
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_1} --includebed ~{truth_bed} -b ~{training_resource_old_08x_vcf} -c ~{training_resource_25x_vcf} -o truvari/
        mv truvari/summary.json ./25x_vs_8x_1.json
        rm -rf truvari/
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_1} --includebed ~{truth_bed} -b ~{training_resource_12x_vcf} -c ~{training_resource_25x_vcf} -o truvari/
        mv truvari/summary.json ./25x_vs_12x_1.json
        rm -rf truvari/
        
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_2} --includebed ~{truth_bed} -b ~{training_resource_old_08x_vcf} -c ~{training_resource_12x_vcf} -o truvari/
        mv truvari/summary.json ./12x_vs_8x_2.json
        rm -rf truvari/
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_2} --includebed ~{truth_bed} -b ~{training_resource_old_08x_vcf} -c ~{training_resource_25x_vcf} -o truvari/
        mv truvari/summary.json ./25x_vs_8x_2.json
        rm -rf truvari/
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS_2} --includebed ~{truth_bed} -b ~{training_resource_12x_vcf} -c ~{training_resource_25x_vcf} -o truvari/
        mv truvari/summary.json ./25x_vs_12x_2.json
        rm -rf truvari/
    >>>

    output {
        Array[File] summaries = glob(work_dir + "/*.json")
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


#
task BenchDipcalls {
    input {
        File dipcall_old_intersample_vcf
        File dipcall_old_intersample_tbi
        File dipcall_asm10_intersample_vcf
        File dipcall_asm10_intersample_tbi
        
        File truth_bed
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 5*ceil(size(dipcall_old_intersample_vcf, "GB") + size(dipcall_asm10_intersample_vcf, "GB")) + 50
    Int ram_gb = 16
    Int min_sv_length = 50
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_RAM_GB=$(( ~{ram_gb} - 2 ))
        TRUVARI_BENCH_SETTINGS="--sizefilt ~{min_sv_length} --sizemin ~{min_sv_length} --sizemax 1000000"
        # Fixed to 50bp since the old dipcall was filtered at 50bp.
        
        ${TIME_COMMAND} truvari bench ${TRUVARI_BENCH_SETTINGS} --includebed ~{truth_bed} -b ~{dipcall_old_intersample_vcf} -c ~{dipcall_asm10_intersample_vcf} -o truvari/
        mv truvari/summary.json ./summary.json
        rm -rf truvari/
    >>>

    output {
        File summary = work_dir + "/summary.json"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 8
        memory: ram_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
