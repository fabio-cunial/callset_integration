version 1.0


# Copies some FORMAT fields of one VCF to corresponding INFO fields of another
# VCF.
#
workflow TransferFormat {
    input {
        File from_vcf_gz
        File from_vcf_gz_tbi
        File to_vcf_gz
        File to_vcf_gz_tbi
        File annotations_map
    }
    parameter_meta {
        annotations_map: "Format: FROM,TO"
    }

    call TransferFormatImpl {
        input:
            from_vcf_gz = from_vcf_gz,
            from_vcf_gz_tbi = from_vcf_gz_tbi,
            to_vcf_gz = to_vcf_gz,
            to_vcf_gz_tbi = to_vcf_gz_tbi,
            annotations_map = annotations_map
    }
    
    output {
        File out_vcf = TransferFormatImpl.out_vcf
        File out_tbi = TransferFormatImpl.out_tbi
    }
}


task TransferFormatImpl {
    input {
        File from_vcf_gz
        File from_vcf_gz_tbi
        File to_vcf_gz
        File to_vcf_gz_tbi
        File annotations_map
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 10*( ceil(size(from_vcf_gz,"GB")) + ceil(size(to_vcf_gz,"GB")) )

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        FORMAT='%CHROM\t%POS\t%ID\t%REF\t%ALT'
        COLUMNS="CHROM,POS,ID,REF,ALT"
        touch header.txt
        while read ROW; do
            ANNOTATION_FROM=${ROW%,*}
            ANNOTATION_TO=${ROW#*,}
            bcftools view --header-only ~{from_vcf_gz} | grep ${ANNOTATION_FROM} | sed -e "s/##FORMAT=<ID=${ANNOTATION_FROM}/##INFO=<ID=${ANNOTATION_TO}/g" >> header.txt
            FORMAT="${FORMAT}\t[%${ANNOTATION_FROM}]"
            COLUMNS="${COLUMNS},INFO/${ANNOTATION_TO}"
        done < ~{annotations_map}
        bcftools query ~{from_vcf_gz} --format ${FORMAT} | bgzip -c > annotations.tsv.gz
        tabix -f -s1 -b2 -e2 annotations.tsv.gz
        bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns ${COLUMNS} --output-type z ~{to_vcf_gz} > out.vcf.gz
        tabix -f out.vcf.gz
    >>>

    output {
        File out_vcf = work_dir + "/out.vcf.gz"
        File out_tbi = work_dir + "/out.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 1
        memory: "8GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}