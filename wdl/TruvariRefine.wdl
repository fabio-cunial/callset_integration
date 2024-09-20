version 1.0


workflow TruvariRefine {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File reference_fa
        File reference_fai
        File? includebed
        String truvari_bench_args = ""
        String truvari_refine_args = ""
    }
    parameter_meta {
        truvari_bench_args: "--sizefilt minSVLengthTest --sizemin minSVLengthTruth"
        truvari_refine_args: "--mafft-params '--auto --thread 8' "
    }
    
    call TruvariRefineImpl {
        input:
            input_vcf_gz = input_vcf_gz,
            input_vcf_gz_tbi = input_vcf_gz_tbi,
            truth_vcf_gz = truth_vcf_gz,
            truth_vcf_gz_tbi = truth_vcf_gz_tbi,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            includebed = includebed,
            truvari_bench_args = truvari_bench_args,
            truvari_refine_args = truvari_refine_args
    }
    
    output {
    	File out = TruvariRefineImpl.out
    }
}


# Remark: the updated TP/FP stats are in $truvari_output/refine.*$.
# $truvari_output/phab_bench$ contains just results from the subset of regions
# which were harmonized.
#
task TruvariRefineImpl {
    input {
        File input_vcf_gz
        File input_vcf_gz_tbi
        File truth_vcf_gz
        File truth_vcf_gz_tbi
        File reference_fa
        File reference_fai
        File? includebed
        String truvari_bench_args = ""
        String truvari_refine_args = ""
    }
    parameter_meta {
    }
    
    String docker_dir = "/truvari_intrasample"
    String work_dir = "/cromwell_root/truvari_intrasample"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        if ~{defined(includebed)} ; then
            INCLUDE_BED="--includebed ~{includebed}"
        else
            INCLUDE_BED=""
        fi
        ${TIME_COMMAND} truvari bench ${INCLUDE_BED} -f ~{reference_fa} -b ~{truth_vcf_gz} -c ~{input_vcf_gz} ~{truvari_bench_args} -o ./truvari_output/
        ${TIME_COMMAND} truvari refine --threads ${N_THREADS} --use-region-coords --use-original-vcfs --recount -f ~{reference_fa} --regions ./truvari_output/candidate.refine.bed ~{truvari_refine_args} ./truvari_output/
        tar -czf out.tar.gz ./truvari_output
    >>>
    
    output {
    	File out = "~{work_dir}/out.tar.gz"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 32
        memory: "128GB"
        disks: "local-disk 100 HDD"
        preemptible: 0
    }
}
