version 1.0


# Remark: the output VCF has a single SAMPLE column that is the concatenation
# of all 5 sample columns emitted by lrcaller (corresponding to different
# genotyping methods, e.g. joint, ad, va, multi).
#
workflow GetRegenotypedVcfLrcaller {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int svlen_min
        Int svlen_max
    }
    parameter_meta {
    }
    
    call GetRegenotypedVcfImpl {
        input:
            truvari_collapsed_vcf_gz = truvari_collapsed_vcf_gz,
            truvari_collapsed_tbi = truvari_collapsed_tbi,
            alignments_bam = alignments_bam,
            alignments_bai = alignments_bai,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            svlen_min = svlen_min,
            svlen_max = svlen_max
    }
    
    output {
        File regenotyped_lrcaller = GetRegenotypedVcfImpl.regenotyped_lrcaller
        File regenotyped_lrcaller_tbi = GetRegenotypedVcfImpl.regenotyped_lrcaller_tbi
    }
}


task GetRegenotypedVcfImpl {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int svlen_min
        Int svlen_max
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    String output_prefix = "lrcaller_regenotyped"
    Int disk_size_gb = 500 + ceil(size(reference_fa,"GB")) + 10*ceil(size(truvari_collapsed_vcf_gz,"GB")) + ceil(size(alignments_bam,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))


        # Makes sure that the merged VCF is in the right format and contains
        # only a specific set of calls.
        #
        function formatVcf() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF_GZ=$2
            
            # Removing multiallelic records
            rm -f tmp0.vcf.gz*
            bcftools norm --multiallelics - --output-type z ${INPUT_VCF_GZ} > tmp0.vcf.gz
            tabix -f tmp0.vcf.gz
            
            # Adding REF and ALT with Adam's script, and removing BNDs.
            rm -f tmp1.vcf*
            python ~{docker_dir}/resolve.py tmp0.vcf.gz ~{reference_fa} | bcftools view -i "SVTYPE!='BND'" | bcftools sort -O z -o tmp1.vcf.gz
            tabix -f tmp1.vcf.gz

            # Keeping only calls in the given length range.
            rm -rf tmp2.vcf*
            FILTER_STRING="SVLEN>=~{svlen_min} && SVLEN<=~{svlen_max}"
            bcftools filter -i "${FILTER_STRING}" --output-type z tmp1.vcf.gz > tmp2.vcf.gz
            tabix -f tmp2.vcf.gz

            # Finalizing
            cp tmp2.vcf.gz ${OUTPUT_VCF_GZ}
            cp tmp2.vcf.gz.tbi ${OUTPUT_VCF_GZ}.tbi
            rm -f tmp*.vcf*
        }
        
        
        # Unpacks truvari's SUPP field into 3 INFO tags.
        #
        function transferSupp() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_VCF_GZ=$2
            
            bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\t[%SUPP]\n' ${INPUT_VCF_GZ} | awk 'BEGIN { FS="\t"; OFS="\t"; } { \
                printf("%s",$1); \
                for (i=2; i<=NF-1; i++) printf("\t%s",$i); \
                if ($6=="0") printf("\t0\t0\t0");
                else if ($6=="1") printf("\t0\t0\t1");
                else if ($6=="2") printf("\t0\t1\t0");
                else if ($6=="3") printf("\t0\t1\t1");
                else if ($6=="4") printf("\t1\t0\t0");
                else if ($6=="5") printf("\t1\t0\t1");
                else if ($6=="6") printf("\t1\t1\t0");
                else if ($6=="7") printf("\t1\t1\t1");
                printf("\n"); \
            }' | bgzip -c > annotations.tsv.gz
            tabix -f -s1 -b2 -e2 annotations.tsv.gz
            rm -f header.txt
            echo '##INFO=<ID=SUPP_PAV,Number=1,Type=Integer,Description="Supported by PAV">' > header.txt
            echo '##INFO=<ID=SUPP_SNIFFLES,Number=1,Type=Integer,Description="Supported by Sniffles2">' >> header.txt
            echo '##INFO=<ID=SUPP_PBSV,Number=1,Type=Integer,Description="Supported by pbsv">' >> header.txt
            bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,INFO/SUPP_PAV,INFO/SUPP_SNIFFLES,INFO/SUPP_PBSV ${INPUT_VCF_GZ} --output-type z > ${OUTPUT_VCF_GZ}
            tabix -f ${OUTPUT_VCF_GZ}
        }


        # Main program
        rm -f ~{alignments_bai}
        samtools index -@ ${N_THREADS} ~{alignments_bam}
    
        # Formatting the merged VCF
        HAS_SUPP=$(bcftools view --header-only ~{truvari_collapsed_vcf_gz} | grep '##FORMAT=<ID=SUPP,' && echo 1 || echo 0)
        if [ ${HAS_SUPP} -eq 0 ]; then
            mv ~{truvari_collapsed_vcf_gz} tmp.vcf.gz
            mv ~{truvari_collapsed_tbi} tmp.vcf.gz.tbi
        else
            transferSupp ~{truvari_collapsed_vcf_gz} tmp.vcf.gz
        fi    
        formatVcf tmp.vcf.gz merged.vcf.gz

        # LRCALLER
        gunzip -c merged.vcf.gz | cut -f 1-8 > merged_cleaned.vcf
        bgzip merged_cleaned.vcf
        tabix -f merged_cleaned.vcf.gz
        ${TIME_COMMAND} lrcaller --number_of_threads ${N_THREADS} --dyn-w-size --fa ~{reference_fa} ~{alignments_bam} merged_cleaned.vcf.gz tmp.vcf 2> /dev/null
        grep '^[#]' tmp.vcf > tmp.txt
        N_ROWS=$(wc -l < tmp.txt)
        head -n $(( ${N_ROWS} - 1 )) tmp.txt > ~{output_prefix}.vcf
        echo '##FORMAT=<ID=GT1,Number=1,Type=String,Description="Genotype">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=GT2,Number=1,Type=String,Description="Genotype">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=GT3,Number=1,Type=String,Description="Genotype">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=GT4,Number=1,Type=String,Description="Genotype">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=GT5,Number=1,Type=String,Description="Genotype">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=AD1,Number=3,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=AD2,Number=3,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=AD3,Number=3,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=AD4,Number=3,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=AD5,Number=3,Type=Integer,Description="Allelic depths from alignment supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=VA1,Number=3,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=VA2,Number=3,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=VA3,Number=3,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=VA4,Number=3,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=VA5,Number=3,Type=Integer,Description="Allelic depths from bam file supporting ref and alt allele and total number of reads">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=PL1,Number=G,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=PL2,Number=G,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=PL3,Number=G,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=PL4,Number=G,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ~{output_prefix}.vcf
        echo '##FORMAT=<ID=PL5,Number=G,Type=Integer,Description="PHRED-scaled genotype likelihoods">' >> ~{output_prefix}.vcf
        echo -e '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE' >> ~{output_prefix}.vcf
        grep '^[^#]' tmp.vcf | awk '{ \
            new_format="GT1:AD1:VA1:PL1:GT2:AD2:VA2:PL2:GT3:AD3:VA3:PL3:GT4:AD4:VA4:PL4:GT5:AD5:VA5:PL5"; \
            new_gt = $10 ":" $11 ":" $12 ":" $13 ":" $14 ; \
            printf("%s",$1); \
            for (i=2; i<=8; i++) printf("\t%s",$i); \
            printf("\t%s",new_format); \
            printf("\t%s",new_gt); \
            printf("\n"); \
        }' >> ~{output_prefix}.vcf
        bgzip ~{output_prefix}.vcf
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_lrcaller = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_lrcaller_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "32GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
