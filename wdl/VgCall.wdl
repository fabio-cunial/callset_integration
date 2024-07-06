version 1.0


# The output by `vg call` contains the following SV features:
#
# INFO.AT="Allele Traversal as path in graph"
# INFO.DP,="Total Depth"
# FORMAT.AD="Allelic depths for the ref and alt alleles in the order listed"
# FORMAT.MAD="Minimum site allele depth"
# FORMAT.DP="Read Depth"
# FORMAT.GL="Genotype Likelihood, log10-scaled likelihoods of the data given
# the called genotype for each possible genotype generated from the reference
# and alternate alleles given the sample ploidy"
# FORMAT.GQ="Genotype Quality, the Phred-scaled probability estimate of the
# called genotype"
# FORMAT.GP="Genotype Probability, the log-scaled posterior probability of the
# called genotype"
# FORMAT.XD="eXpected Depth, background coverage as used for the Poisson model"
# FILTER.lowad="Variant does not meet minimum allele read support threshold of
# 1"
# FILTER.lowdepth="Variant has read depth less than 4"
#
workflow VgCall {
    input {
        String remote_dir
        String sample_id
        File alignments_gam
        File vcf_gz
        File vcf_tbi
        Int min_mapq = 5
        Int trim_ends = 5
        Int strict = 1
        Int n_cpus = 32
        Int ram_size_gb = 128
        Int disk_size_gb = 256
    }
    parameter_meta {
        remote_dir: "The XG graph index is assumed to be in directory `remote_dir/sample_id`."
        strict: "1=Re-genotypes the same VCF that was used to build the graph. 0=Uses the graph built from the original VCF to call SVs. This may output a different VCF."
    }
    call VgCallImpl {
        input:
            remote_dir = remote_dir,
            sample_id = sample_id,
            alignments_gam = alignments_gam,
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            min_mapq = min_mapq,
            trim_ends = trim_ends,
            strict = strict,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File genotyped_vcf = VgCallImpl.genotyped_vcf
        File genotyped_tbi = VgCallImpl.genotyped_tbi
    }
}


# COMMAND        | TIME | CORES | RAM
# vg pack        | 1h   |  27   | 61 G
# vg call --vcf  | 1h   |  21   | 128 G
# vg call        | 1h   |  30   | 115 G
#
task VgCallImpl {
    input {
        String remote_dir
        String sample_id
        File alignments_gam
        File vcf_gz
        File vcf_tbi
        Int min_mapq
        Int trim_ends
        Int strict
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
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
        VG_COMMAND="~{docker_dir}/vg"
        
        if [ ~{strict} -eq 1 ]; then
            VCF_FLAG="--vcf ~{vcf_gz}"
        else
            VCF_FLAG=""
        fi
        while : ; do
            TEST=$(gsutil -m cp ~{remote_dir}/~{sample_id}/~{sample_id}.xg . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading XG. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        touch ~{vcf_tbi}
        ${TIME_COMMAND} ${VG_COMMAND} pack --threads ${N_THREADS} --xg ~{sample_id}.xg --gam ~{alignments_gam} --min-mapq ~{min_mapq} --trim-ends ~{trim_ends} --packs-out ~{sample_id}.pack
        ${TIME_COMMAND} ${VG_COMMAND} call --threads ${N_THREADS} --progress --pack ~{sample_id}.pack --ploidy 2 ${VCF_FLAG} ~{sample_id}.xg > out.vcf
        df -h
        bcftools sort -m $(( ~{ram_size_gb} - 8 ))G --output-type z > ~{sample_id}.vcf.gz out.vcf
        tabix -f ~{sample_id}.vcf.gz
        df -h
    >>>
    
    output {
        File genotyped_vcf = work_dir + "/" + sample_id + ".vcf.gz"
        File genotyped_tbi = work_dir + "/" + sample_id + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
