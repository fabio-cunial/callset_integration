version 1.0


#
workflow GiraffeMap {
    input {
        String remote_dir
        String sample_id
        File reads1_fastq_gz
        File reads2_fastq_gz
        Int n_cpus = 16
        Int ram_size_gb = 32
        Int disk_size_gb
    }
    parameter_meta {
        remote_dir: "The GBZ, MIN and DIST graph indexes are assumed to be in directory `remote_dir/sample_id`."
    }
    call GiraffeMapImpl {
        input:
            remote_dir = remote_dir,
            sample_id = sample_id,
            reads1_fastq_gz = reads1_fastq_gz,
            reads2_fastq_gz = reads2_fastq_gz,
            n_cpus = n_cpus,
            ram_size_gb = ram_size_gb,
            disk_size_gb = disk_size_gb
    }
    output {
        File alignments_gam = GiraffeMapImpl.alignments_gam
        File stats = GiraffeMapImpl.stats
    }
}


# COMMAND    | TIME | CORES | RAM
# vg giraffe | 6h   |   16  | 50G
# vg stats   |  ?   |    ?  |   ?     CRASHED
#
task GiraffeMapImpl {
    input {
        String remote_dir
        String sample_id
        File reads1_fastq_gz
        File reads2_fastq_gz
        Int n_cpus
        Int ram_size_gb
        Int disk_size_gb
    }
    parameter_meta {
    }
    
    String docker_dir = "/infogain"
    String work_dir = "/cromwell_root/infogain"
    
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
        
        
        while : ; do
            TEST=$(gsutil -m cp ~{remote_dir}/~{sample_id}/~{sample_id}.gbz ~{remote_dir}/~{sample_id}/~{sample_id}.min ~{remote_dir}/~{sample_id}/~{sample_id}.dist . && echo 0 || echo 1)
            if [ ${TEST} -eq 1 ]; then
                echo "Error downloading indexes. Trying again..."
                sleep ${GSUTIL_DELAY_S}
            else
                break
            fi
        done
        ${TIME_COMMAND} ${VG_COMMAND} giraffe --threads ${N_THREADS} --progress --output-format gam --gbz-name ~{sample_id}.gbz --minimizer-name ~{sample_id}.min --dist-name ~{sample_id}.dist --fastq-in ~{reads1_fastq_gz} --fastq-in ~{reads2_fastq_gz} > ~{sample_id}.gam
        # The following command crashes on vg 1.58.0 (full error log below):
        #${TIME_COMMAND} ${VG_COMMAND} stats --threads ${N_THREADS} --alignments ~{sample_id}.gam > stats.txt
        touch stats.txt
        df -h
    >>>
    
    output {
        File alignments_gam = work_dir + "/" + sample_id + ".gam"
        File stats = work_dir + "/stats.txt"
    }
    runtime {
        docker: "fcunial/infogain"
        cpu: n_cpus
        memory: ram_size_gb + "GB"
        disks: "local-disk " + disk_size_gb + " HDD"
        preemptible: 0
    }
}


# vg stats error:
#
# terminate called after throwing an instance of 'std::runtime_error'
#   what():  [vg::io::MessageIterator] obsolete, invalid, or corrupt input at message 9202333682039739 group 9202327764364411
# ━━━━━━━━━━━━━━━━━━━━
# Crash report for vg v1.58.0 "Cartari"
# Stack trace (most recent call last):
# #16   Object "/infogain/vg", at 0x61f764, in _start
# #15   Object "/infogain/vg", at 0x20f7be6, in __libc_start_main
# #14   Object "/infogain/vg", at 0x20f6349, in __libc_start_call_main
# #13   Object "/infogain/vg", at 0xe14beb, in vg::subcommand::Subcommand::operator()(int, char**) const
# #12   Object "/infogain/vg", at 0xe0dd92, in main_stats(int, char**)
# #11   Object "/infogain/vg", at 0x20d3b05, in GOMP_parallel
# #10   Object "/infogain/vg", at 0xe0b87c, in void vg::io::for_each_parallel_impl<vg::Alignment>(std::istream&, std::function<void (vg::Alignment&, vg::Alignment&)> const&, std::function<void (vg::Alignment&)> const&, std::function<bool ()> const&, unsigned long) [clone ._omp_fn.0]
# #9    Object "/infogain/vg", at 0x5dd7d1, in vg::io::MessageIterator::take[abi:cxx11]() [clone .cold]
# #8    Object "/infogain/vg", at 0x20f3a4d, in _Unwind_Resume
# #7    Object "/infogain/vg", at 0x20f2ecb, in _Unwind_RaiseException_Phase2
# #6    Object "/infogain/vg", at 0x202f689, in __gxx_personality_v0
# #5    Object "/infogain/vg", at 0x20c9b98, in __cxa_call_terminate
# #4    Object "/infogain/vg", at 0x202ff3b, in __cxxabiv1::__terminate(void (*)())
# #3    Object "/infogain/vg", at 0x5eba23, in __gnu_cxx::__verbose_terminate_handler() [clone .cold]
# #2    Object "/infogain/vg", at 0x5ee16b, in abort
# #1    Object "/infogain/vg", at 0x210edc5, in raise
# #0    Object "/infogain/vg", at 0x213b91c, in __pthread_kill
# ERROR: Signal 6 occurred. VG has crashed. Visit https://github.com/vgteam/vg/issues/new/choose to report a bug.
# Please include this entire error log in your bug report!
# ━━━━━━━━━━━━━━━━━━━━
# Command exited with non-zero status 134
#     Command being timed: "/infogain/vg stats --threads 16 --alignments HG002.gam"
#     User time (seconds): 45070.25
#     System time (seconds): 438.51
#     Percent of CPU this job got: 1449%
#     Elapsed (wall clock) time (h:mm:ss or m:ss): 52:19.45
#     Average shared text size (kbytes): 0
#     Average unshared data size (kbytes): 0
#     Average stack size (kbytes): 0
#     Average total size (kbytes): 0
#     Maximum resident set size (kbytes): 61021544
#     Average resident set size (kbytes): 0
#     Major (requiring I/O) page faults: 154
#     Minor (reclaiming a frame) page faults: 22166815
#     Voluntary context switches: 12789848
#     Involuntary context switches: 7646085
#     Swaps: 0
#     File system inputs: 274279736
#     File system outputs: 8
#     Socket messages sent: 0
#     Socket messages received: 0
#     Signals delivered: 0
#     Page size (bytes): 4096
#     Exit status: 134