version 1.0


#
workflow AouVsHPRC {
    input {
        File aou_vcf_gz
        File hprc_dipcall_vcf_gz
        String truvari_args = "--sizefilt 0 --sizemax 1000000 --sizemin 0 --pctsize 0.9 --pctseq 0"
    }
    parameter_meta {
    }

    call AouVsHPRCImpl {
        input:
            aou_vcf_gz = aou_vcf_gz,
            hprc_dipcall_vcf_gz = hprc_dipcall_vcf_gz,
            truvari_args = truvari_args
    }
    
    output {
        File truvari_output = AouVsHPRCImpl.truvari_output
    }
}


task AouVsHPRCImpl {
    input {
        File aou_vcf_gz
        File hprc_dipcall_vcf_gz
        String truvari_args
    }
    parameter_meta {
    }
    
    String work_dir = "/cromwell_root/aou"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( 2 * ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        
        
        function format() {
            local INPUT_VCF_GZ=$1
            local OUTPUT_FILE=$2
            
            ${TIME_COMMAND} bcftools view --threads ${N_THREADS} --drop-genotypes --output-type z ${INPUT_VCF_GZ > tmp.vcf.gz
            tabix tmp.vcf.gz
            bcftools view --header-only tmp.vcf.gz > header.txt
            N_ROWS=$(wc -l < header.txt)
            head -n $(( ${N_ROWS} - 1 )) header.txt > out.vcf
            echo "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" >> out.vcf
            bcftools view --no-header tmp.vcf.gz | awk '{ printf("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tGT\t0/1\n",$1,$2,$3,$4,$5,$6,$7,$8); }' >> out.vcf
            bgzip -c out.vcf > ${OUTPUT_FILE}
            tabix ${OUTPUT_FILE}
            rm -f tmp.vcf* header.txt out.vcf*
        }
        
        
        format ~{aou_vcf_gz} aou.vcf.gz
        format ~{hprc_dipcall_vcf_gz} hprc.vcf.gz
        ${TIME_COMMAND} truvari bench ~{truvari_args} -b hprc.vcf.gz -c aou.vcf.gz -o ./truvari_dir
        tar -czf output.tar.gz ./truvari_dir/
    >>>
    
    output {
        File truvari_output = work_dir + "/output.tar.gz"
    }
    runtime {
        docker: "us.gcr.io/broad-dsp-lrma/aou-lr/truvari_intrasample"
        cpu: 8
        memory: "32GB"
        disks: "local-disk 256 HDD"
        preemptible: 0
    }
}