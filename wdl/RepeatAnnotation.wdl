version 1.0


# Adds a new field INFO/tr_coverage that contains the fraction of a call covered
# by intervals in a given BED.
#
workflow RepeatAnnotation {
    input {
        File vcf_gz_file
        File tbi_file
        File bed_file
    }
    parameter_meta {
    }

    call RepeatAnnotationImpl {
        input:
            vcf_gz_file = vcf_gz_file,
            tbi_file = tbi_file,
            bed_file = bed_file
    }
    
    output {
        File vcf_gz = RepeatAnnotationImpl.annotated_vcf
        File vcf_gz_tbi = RepeatAnnotationImpl.annotated_tbi
    }
}


task RepeatAnnotationImpl {
    input {
        File vcf_gz_file
        File tbi_file
        File bed_file
    }
    parameter_meta {
    }
    
    String docker_dir = "/callset_integration"
    String work_dir = "/cromwell_root/callset_integration"
    Int disk_size_gb = 10*ceil(size(vcf_gz_file,"GB")) + ceil(size(bed_file,"GB"))

    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        echo '##INFO=<ID=tr_coverage,Number=1,Type=Float,Description="Fraction of a call that is covered by intervals in a TR BED.">' > header.txt
        bcftools query --format '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' ~{vcf_gz_file} > tmp1.tsv
        bedtools annotate -i ~{vcf_gz_file} -files ~{bed_file} | cut -f 1,2,3,4,5,11 | sort -k 1V -k 2n | bgzip -c > annotations.tsv.gz
        tabix -f -s1 -b2 -e2 annotations.tsv.gz
        bcftools annotate --annotations annotations.tsv.gz --header-lines header.txt --columns CHROM,POS,ID,REF,ALT,INFO/tr_coverage ~{vcf_gz_file} --output-type z > out.vcf.gz
        tabix -f out.vcf.gz
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