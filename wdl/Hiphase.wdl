version 1.0


# 
#
workflow Hiphase {
    input {
        String samplename
        File unphased_vcf
        File? regions_bed
        File bam
        File bai
        File ref_fa
        File ref_fai
        Int max_sv_length = 99000
    }
    parameter_meta {
        regions_bed: "Subset the VCF to these regions."
        max_sv_length: "To avoid the following error: 'Encountered WFA error for mapping [...]: Max_edit_distance (100000) reached during WFA solving'."
    }

    call HiphaseImpl {
        input:
            samplename = samplename,
            unphased_vcf = unphased_vcf,
            regions_bed = regions_bed,
            bam = bam,
            bai = bai,
            ref_fa = ref_fa,
            ref_fai = ref_fai,
            max_sv_length = max_sv_length
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
        File? regions_bed
        File bam
        File bai
        File ref_fa
        File ref_fai
        Int max_sv_length
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
        
        # Keeping only calls in the BED file (if any). 
        tabix -f ~{unphased_vcf}
        if ~{defined(regions_bed)} ; then
            bcftools view --regions-file ~{regions_bed} --output-type z ~{unphased_vcf} > tmp1.vcf.gz
            tabix -f tmp1.vcf.gz
        else
            mv ~{unphased_vcf} tmp1.vcf.gz
            mv ~{unphased_vcf}.tbi tmp1.vcf.gz.tbi
        fi
            
        # Ensuring the right sample name in the VCF.
        echo ~{samplename} > sample.txt
        bcftools reheader --threads ${N_THREADS} --samples sample.txt --output tmp2.vcf.gz tmp1.vcf.gz
        tabix -f tmp2.vcf.gz
        
        # Ensuring that the reference characters agree with the FASTA.
        bcftools norm --threads ${N_THREADS} --check-ref s --fasta-ref ~{ref_fa} --do-not-normalize --output-type z tmp2.vcf.gz > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        
        # Ensuring all uppercases in REF, ALT, and reference FASTA.
        # Keeping only calls up to the given max length.
        bcftools view --header-only tmp3.vcf.gz > filtered.vcf
        bcftools view --no-header tmp3.vcf.gz | awk '{ \
            if (length($4)<=~{max_sv_length} && length($5)<=~{max_sv_length}) { \
                $4=toupper($4); \
                if (substr($5,1,1)!="<") $5=toupper($5); \
                printf("%s",$1); \
                for (i=2; i<=NF; i++) printf("\t%s",$i); \
                printf("\n"); \
            } \
        }' >> filtered.vcf
        bgzip filtered.vcf; tabix -f filtered.vcf.gz
        bcftools view --no-header filtered.vcf.gz | head || echo 1
        awk '{ \
            if (substr($0,1,1)!=">") $0=toupper($0); \
            printf("%s",$0);
            printf("\n"); \
        }' ~{ref_fa} > reference.fa
        rm -f ~{ref_fa}
        mv ~{ref_fai} reference.fa.fai
        head reference.fa || echo 1
        
        # Phasing
        touch ~{bai}
        ${TIME_COMMAND} hiphase \
            --threads ${N_THREADS} \
            --bam ~{bam} \
            --reference reference.fa \
            --global-realignment-cputime 300 \
            --vcf filtered.vcf.gz \
            --output-vcf ~{samplename}_phased.vcf.gz \
            --output-bam ~{samplename}_haplotagged.bam \
            --haplotag-file ~{samplename}.haplotag.tsv \
            --stats-file ~{samplename}.stats.csv \
            --blocks-file ~{samplename}.blocks.tsv \
            --summary-file ~{samplename}.summary.tsv
        ${TIME_COMMAND} bcftools sort ~{samplename}_phased.vcf.gz -O z -o ~{samplename}_phased.sorted.vcf.gz
        tabix -f ~{samplename}_phased.sorted.vcf.gz
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
