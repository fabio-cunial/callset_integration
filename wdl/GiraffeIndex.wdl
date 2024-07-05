version 1.0


#
workflow GiraffeIndex {
    input {
        String sample_id
        File input_vcf_gz
        File input_tbi
        String remote_dir
        File reference_fa
        File reference_fai
        Int n_cpus = 16
        Int ram_size_gb = 128
        Int disk_size_gb = 512
    }
    parameter_meta {
        remote_dir: "The program stores all graph indexes in `remote_dir/sample_id`."
    }
    call GiraffeIndexImpl {
        input:
            sample_id = sample_id,
            input_vcf_gz = input_vcf_gz,
            input_tbi = input_tbi,
            remote_dir = remote_dir,
            reference_fa = reference_fa,
            reference_fai = reference_fai,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    output {
    }
}


# See this GitHub issue <https://github.com/vgteam/vg/issues/3950> for details
# about how to build the indexes and on how to run `vg call`.
#
# COMMAND              | TIME  | CORES | RAM
# vg construct         | 6m    | 1     | 54 M
# vg index --dist-name | 1h30m | 16    | 104 G
# vg index --xg-name   | 1h15m | 1     | 33.5 G
# vg gbwt              | 4h    | 12    | 50 G
# vg gbwt --gbz-format | 3m    | 1     | 18 G
# vg minimizer         | 10m   | 19    | 64 G
#
task GiraffeIndexImpl {
    input {
        String sample_id
        File input_vcf_gz
        File input_tbi
        String remote_dir
        File reference_fa
        File reference_fai
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    Int ram_size_gb_vg = ram_size_gb - 8
    
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
        
        ${TIME_COMMAND} ${VG_COMMAND} construct --threads ${N_THREADS} --progress --handle-sv --alt-paths --reference ~{reference_fa} --vcf ~{input_vcf_gz} > ~{sample_id}.vg
        mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} index --threads ${N_THREADS} --temp-dir ./vgtmp --progress --dist-name ~{sample_id}.dist ~{sample_id}.vg
        rm -rf ./vgtmp; mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} index --threads ${N_THREADS} --temp-dir ./vgtmp --progress --xg-alts --xg-name ~{sample_id}.xg ~{sample_id}.vg
        rm -rf ./vgtmp; mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} gbwt --num-jobs ${N_THREADS} --temp-dir ./vgtmp --progress --path-cover --xg-name ~{sample_id}.xg --output ~{sample_id}.gbwt
        rm -rf ./vgtmp; mkdir ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} gbwt --num-jobs ${N_THREADS} --temp-dir ./vgtmp --progress --xg-name ~{sample_id}.xg --graph-name ~{sample_id}.gbz --gbz-format ~{sample_id}.gbwt
        rm -rf ./vgtmp
        ${TIME_COMMAND} ${VG_COMMAND} minimizer --threads ${N_THREADS} --progress --distance-index ~{sample_id}.dist --output-name ~{sample_id}.min ~{sample_id}.gbz
        df -h
        
        # Uploading all indexes
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.xg ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading XG. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.gbz ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading GBZ. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
         while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.min ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading MIN. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.dist ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading DIST. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
         while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.vg ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading VG. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        while : ; do
            TEST=$(gsutil ${GSUTIL_UPLOAD_THRESHOLD} -m cp ~{sample_id}.gbwt ~{remote_dir}/~{sample_id}/ && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error uploading GBWT. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
    >>>
    
    output {
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}
