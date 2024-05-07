version 1.0


#
workflow GetRegenotypedVcfCutesv {
    input {
        File truvari_collapsed_vcf_gz
        File truvari_collapsed_tbi
        File alignments_bam
        File alignments_bai
        File reference_fa
        File reference_fai
        Int svlen_min
        Int svlen_max
        String cutesv_params = "--max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 --merge_ins_threshold 500 --merge_del_threshold 500"
        Int cutesv_minsupport = 10
    }
    parameter_meta {
        cutesv_params: "Default values are taken from the command-line help."
        cutesv_minsupport: "Set it to one if the purpose is just annotating every call with features. Default=10."
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
            svlen_max = svlen_max,
            cutesv_params = cutesv_params,
            cutesv_minsupport = cutesv_minsupport
    }
    
    output {
        File regenotyped_cutesv = GetRegenotypedVcfImpl.regenotyped_cutesv
        File regenotyped_cutesv_tbi = GetRegenotypedVcfImpl.regenotyped_cutesv_tbi
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
        String cutesv_params
        Int cutesv_minsupport
    }
    parameter_meta {
    }
    
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    String output_prefix = "cutesv_regenotyped"
    Int disk_size_gb = 250 + ceil(size(reference_fa,"GB")) + 10*ceil(size(truvari_collapsed_vcf_gz,"GB")) + ceil(size(alignments_bam,"GB"))

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

        # CUTESV FORCE
        rm -rf ./cutesv_tmp; mkdir ./cutesv_tmp
        ${TIME_COMMAND} cuteSV --threads ${N_THREADS} -Ivcf merged.vcf.gz ~{cutesv_params} --genotype --min_support ~{cutesv_minsupport} --min_size ~{svlen_min} --max_size ~{svlen_max} ~{alignments_bam} ~{reference_fa} tmp1.vcf ./cutesv_tmp
        rm -rf ./cutesv_tmp
        bgzip tmp1.vcf
        tabix -f tmp1.vcf.gz
        
        # Fixing SAMPLE (set to NULL by cutesv)
        bcftools view --header-only merged.vcf.gz | tail -n 1 | cut -f 10 > sample.txt
        bcftools reheader --threads ${N_THREADS} --samples sample.txt tmp1.vcf.gz > tmp2.vcf.gz
        tabix -f tmp2.vcf.gz
        rm -f tmp1.vcf.gz*
        
        # Importing cutesv's FORMAT into the original file. This is needed since
        # cutesv destroys several INFO field annotations.
        echo -e "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"# Phred-scaled genotype likelihoods rounded to the closest integer\">" > header.txt
        # The other FORMAT fields are assumed to be already in the input VCF
        bcftools annotate --header-lines header.txt --annotations tmp2.vcf.gz --columns CHROM,POS,ID,REF,ALT,FORMAT/GT,FORMAT/DR,FORMAT/DV,FORMAT/PL,FORMAT/GQ --output-type z merged.vcf.gz > ~{output_prefix}.vcf.gz    
        tabix -f ~{output_prefix}.vcf.gz
    >>>

    output {
        File regenotyped_cutesv = work_dir + "/" + output_prefix + ".vcf.gz"
        File regenotyped_cutesv_tbi = work_dir + "/" + output_prefix + ".vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 16
        memory: "32GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
