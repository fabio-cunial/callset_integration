version 1.0


# 
#
workflow Hiphase {
    input {
        String samplename
        File unphased_vcf
        File unphased_tbi
        File? regions_bed
        File bam
        File bai
        File ref_fa
        File ref_fai
    }
    parameter_meta {
        regions_bed: "Subset the VCF to these regions."
    }

    call HiphaseImpl {
        input:
            samplename = samplename,
            unphased_vcf = unphased_vcf,
            unphased_tbi = unphased_tbi,
            regions_bed = regions_bed,
            bam = bam,
            bai = bai,
            ref_fa = ref_fa,
            ref_fai = ref_fai
    }
    
    output {
        File phased_vcf = HiphaseImpl.phased_vcf
        File phased_tbi = HiphaseImpl.phased_tbi
        File haplotag_file = HiphaseImpl.haplotag_file
        File stats_file = HiphaseImpl.stats_file
        File blocks_file = HiphaseImpl.blocks_file
        File summary_file = HiphaseImpl.summary_file
        File haplotagged_bam = HiphaseImpl.haplotagged_bam
    }
}


task HiphaseImpl {
    input {
        String samplename
        File unphased_vcf
        File unphased_tbi
        File? regions_bed
        File bam
        File bai
        File ref_fa
        File ref_fai
    }
    
    Int disk_size_gb = 5*ceil( size(unphased_vcf,"GB") + size(bam,"GB") + size(ref_fa,"GB") )
    String docker_dir = "/"
    String work_dir = "/cromwell_root/hiphase"

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
    
        TIME_COMMAND=" "
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        if ~{defined(regions_bed)} ; then
            bcftools view --regions-file ~{regions_bed} --output-type z ~{unphased_vcf} > filtered.vcf.gz
            tabix filtered.vcf.gz
        else
            mv ~{unphased_vcf} filtered.vcf.gz
            mv ~{unphased_tbi} filtered.vcf.gz.tbi
        fi    
        touch ~{bai}
        ${TIME_COMMAND} hiphase \
            --threads ${N_THREADS} \
            --bam ~{bam} \
            --reference ~{ref_fa} \
            --global-realignment-cputime 300 \
            --vcf filtered.vcf.gz \
            --output-vcf ~{samplename}_phased.vcf.gz \
            --output-bam ~{samplename}_haplotagged.bam \
            --haplotag-file ~{samplename}.haplotag.tsv \
            --stats-file ~{samplename}.stats.csv \
            --blocks-file ~{samplename}.blocks.tsv \
            --summary-file ~{samplename}.summary.tsv \
            --verbose
        ${TIME_COMMAND} bcftools sort ~{samplename}_phased.vcf.gz -O z -o ~{samplename}_phased.sorted.vcf.gz
        tabix ~{samplename}_phased.sorted.vcf.gz
    >>>

    output {
        File phased_vcf = work_dir + "/" + samplename + "_phased.sorted.vcf.gz"
        File phased_tbi = work_dir + "/" + samplename + "_phased.sorted.vcf.gz.tbi"
        File haplotag_file = work_dir + "/" + samplename + ".haplotag.tsv"
        File stats_file = work_dir + "/" + samplename + ".stats.csv"
        File blocks_file = work_dir + "/" + samplename + ".blocks.tsv"
        File summary_file = work_dir + "/" + samplename + ".summary.tsv"
        File haplotagged_bam = work_dir + "/" + samplename + "_haplotagged.bam"
    }
    runtime {
        docker: "hangsuunc/hiphase:1.3.0"
        cpu: 16
        memory: "32GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
