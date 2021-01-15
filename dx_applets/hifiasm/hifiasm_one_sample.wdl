version 1.0

task hifiasm{
    input {
        Array[File] raw_reads
        String? optional_argument
        String output_prefix
    }
    command <<<
        set -e -x -o pipefail

        cmd=(/hifiasm-0.7/hifiasm -o ~{output_prefix} -t $(nproc))
        for file in "~{sep=' ' raw_reads}"
        do
            cmd+=(${file})
        done
        "${cmd[@]}"
        awk '$1 ~/S/ {print ">"$2"\n"$3}' ~{output_prefix}.p_ctg.gfa > ~{output_prefix}.fasta
    >>>
    output {
        File haplotype_resolve_raw_unitig = "${output_prefix}.r_utg.gfa"
        File haplotype_resolve_raw_unitig_no_bubble = "${output_prefix}.p_utg.gfa"
        File primary_assembly_contig_graph = "${output_prefix}.p_ctg.gfa"
        File alternate_assembly_contig_graph = "${output_prefix}.a_ctg.gfa"
        File primary_assembly = "${output_prefix}.fasta"
    }
    runtime {
        docker: "quay.io/chai/hifiasm:0.7.0"
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
