version 1.0


# Adds a new field INFO/tr_coverage that contains the fraction of a call covered
# by intervals in a given BED.
#
workflow RepeatAnnotation {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File tr_bed
        File segdup_bed
        File telomere_bed
        File centromere_bed
    }
    parameter_meta {
    }

    call RepeatAnnotationImpl {
        input:
            vcf_gz = vcf_gz,
            vcf_gz_tbi = vcf_gz_tbi,
            tr_bed = tr_bed,
            segdup_bed = segdup_bed,
            telomere_bed = telomere_bed,
            centromere_bed = centromere_bed
    }
    
    output {
        File annotated_vcf = RepeatAnnotationImpl.annotated_vcf
        File annotated_tbi = RepeatAnnotationImpl.annotated_tbi
    }
}


task RepeatAnnotationImpl {
    input {
        File vcf_gz
        File vcf_gz_tbi
        File tr_bed
        File segdup_bed
        File telomere_bed
        File centromere_bed
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 10*ceil(size(vcf_gz,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        function annotate() {
            local VCF_GZ_FILE=$1
            local BED_FILE=$2
            local INFO_ID=$3
            local DESCRIPTION=$4
            local OUT_FILE=$5
            
            echo -e "##INFO=<ID=${INFO_ID},Number=1,Type=Float,Description=\"Fraction of a call that is covered by intervals in a ${DESCRIPTION} BED.\">" > header.txt
            bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ${VCF_GZ_FILE} > tmp.tsv
            bedtools annotate -i ${VCF_GZ_FILE} -files ${BED_FILE} | cut -f 1,2,3,4,5,11 | sort -k 1V -k 2n | bgzip -c > annotations.tsv.gz
            tabix -f -s1 -b2 -e2 annotations.tsv.gz
            bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,INFO/${INFO_ID} ${VCF_GZ_FILE} --output-type z > ${OUT_FILE}
            tabix -f ${OUT_FILE}
        }
        
        
        # Main program
        annotate ~{vcf_gz} ~{tr_bed} tr_coverage TR out1.vcf.gz
        annotate out1.vcf.gz ~{segdup_bed} segdup_coverage segdup out2.vcf.gz
        rm -f out1.vcf.gz*
        annotate out2.vcf.gz ~{telomere_bed} telomere_coverage telomere out1.vcf.gz
        rm -f out2.vcf.gz*
        annotate out1.vcf.gz ~{centromere_bed} centromere_coverage centromere out.vcf.gz
    >>>

    output {
        File annotated_vcf = work_dir + "/out.vcf.gz"
        File annotated_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}