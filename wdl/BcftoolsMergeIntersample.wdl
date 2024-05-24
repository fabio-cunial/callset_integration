version 1.0


# Used for building a simple cohort VCF to be filtered with the reads of each
# sample, in the `sv_filtering_manuscript` workspace.
#
workflow BcftoolsMergeIntersample {
    input {
        Array[File] vcf_gz
        Array[File] vcf_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    call BcftoolsMergeIntersampleImpl {
        input:
            vcf_gz = vcf_gz,
            vcf_tbi = vcf_tbi,
            reference_fa = reference_fa
    }
    
    output {
    	File bcftools_merged_vcf = BcftoolsMergeIntersampleImpl.bcftools_merged_vcf
    	File bcftools_merged_tbi = BcftoolsMergeIntersampleImpl.bcftools_merged_tbi
    }
}


#
task BcftoolsMergeIntersampleImpl {
    input {
        Array[File] vcf_gz
        Array[File] vcf_tbi
        File reference_fa
    }
    parameter_meta {
    }
    
    Int disk_size_gb = 10*ceil(size(vcf_gz,"GB")) + 50
    Int mem_gb = 64
    String docker_dir = "/hgsvc2"
    String work_dir = "/cromwell_root/hgsvc2"
    
    command <<<
        set -euxo pipefail
        mkdir -p ~{work_dir}
        cd ~{work_dir}
        
        TIME_COMMAND="/usr/bin/time --verbose"
        N_SOCKETS="$(lscpu | grep '^Socket(s):' | awk '{print $NF}')"
        N_CORES_PER_SOCKET="$(lscpu | grep '^Core(s) per socket:' | awk '{print $NF}')"
        N_THREADS=$(( ${N_SOCKETS} * ${N_CORES_PER_SOCKET} ))
        EFFECTIVE_MEM_GB=~{mem_gb}
        EFFECTIVE_MEM_GB=$(( ${EFFECTIVE_MEM_GB} - 4 ))
        
        INPUT_FILES=~{sep=',' vcf_gz}
        echo ${INPUT_FILES} | tr ',' '\n' > list1.txt
        
        # Formatting input files
        mkdir -p ./preprocessed/
        while read FILE; do
            # Removing multiallelic records, if any.
            bcftools norm --multiallelics - --output-type z ${FILE} > tmp1.vcf.gz
            tabix tmp1.vcf.gz
            # Fixing REF=N (caused e.g. by sniffles).
            bcftools norm --check-ref s --fasta-ref ~{reference_fa} --do-not-normalize --output-type z tmp1.vcf.gz > tmp2.vcf.gz
            tabix -f tmp2.vcf.gz; rm -f tmp1.vcf*
            # Enforcing a consistent format
            outname=./preprocessed/$(basename ${FILE})
            python ~{docker_dir}/resolve_light.py tmp2.vcf.gz ~{reference_fa} | bcftools sort --max-mem ${EFFECTIVE_MEM_GB}G --output-type z > ${outname}
            tabix -f ${outname}; rm -f tmp2.vcf*
            echo ${outname} >> list2.txt
        done < list1.txt

        # Merging
        ${TIME_COMMAND} bcftools merge --threads ${N_THREADS} --merge none --force-samples --file-list list2.txt --output-type z > tmp3.vcf.gz
        tabix -f tmp3.vcf.gz
        
        # Removing multiallelic records again, just to be completely sure they
        # don't exist.
        bcftools norm --multiallelics - --output-type z tmp3.vcf.gz > bcftools_merged.vcf.gz
        tabix -f bcftools_merged.vcf.gz
    >>>
    
    output {
    	File bcftools_merged_vcf = "~{work_dir}/bcftools_merged.vcf.gz"
    	File bcftools_merged_tbi = "~{work_dir}/bcftools_merged.vcf.gz.tbi"
    }
    runtime {
        docker: "fcunial/callset_integration"
        cpu: 2
        memory: mem_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
