version 1.0

task hifiasm{
    input {
        Array[File] raw_reads
        String? optional_argument
        String output_prefix
    }
    command <<<
        set -e -x -o pipefail

        cmd=(/hifiasm-0.15.2/hifiasm -o ~{output_prefix} -t $(nproc))
        for file in "~{sep=' ' raw_reads}"
        do
            cmd+=(${file})
        done
        "${cmd[@]}"
        ls -lhtr
        awk '$1 ~/S/ {print ">"$2"\n"$3}' ~{output_prefix}.bp.p_ctg.gfa > ~{output_prefix}.fasta
        awk '$1 ~/S/ {print ">"$2"\n"$3}' ~{output_prefix}.bp.hap1.p_ctg.gfa > ~{output_prefix}.hap1.fasta
        awk '$1 ~/S/ {print ">"$2"\n"$3}' ~{output_prefix}.bp.hap2.p_ctg.gfa > ~{output_prefix}.hap2.fasta
    >>>
    output {
        File primary_assembly = "${output_prefix}.fasta"
        File hap1 = "${output_prefix}.hap1.fasta"
        File hap2 = "${output_prefix}.hap2.fasta"
        Array[File] all_bed = glob("*.bed")
        Array[File] all_gfa = glob("*.gfa")

    }
    runtime {
        docker: "quay.io/chai/hifiasm:0.15.2"
        cpu: 64
        memory: "500 GB"
        disk: "/tmp 1000 SSD"
        dx_timeout: "120H"
    }
    parameter_meta {
    raw_reads: {
        description: "raw long read data",
        patterns: ["*.fasta","*.fa","*.fastq","*.fq","*.fasta.gz","*.fa.gz","*.fastq.gz","*.fq.gz"],
        stream: true,
        localization_optional: true
    }
    }
}
