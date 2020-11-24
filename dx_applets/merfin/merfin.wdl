version 1.0

task merfin{
    input {
        File readdb_meryl
        File seq_meryl
        File sequence_fasta
        File vcf
        Int peak
    }
    command <<<
        set -e -x -o pipefail
        mkdir readdbmeryl
        mkdir seqdbmeryl
        tar -xvf ~{readdb_meryl} -C readdbmeryl
        tar -xvf ~{seq_meryl} -C seqdbmeryl
        zcat ~{sequence_fasta} > sequence.fasta
        zcat ~{vcf} > vcf.vcf
        mv readdbmeryl/all_count/union_meryl readdbmeryl/all_count/union.meryl
        mv seqdbmeryl/all_count/union_meryl seqdbmeryl/all_count/union.meryl
        /merfin/build/bin/merfin -vmer -memory1 2 -memory2 16 -sequence sequence.fasta -seqmers seqdbmeryl/all_count/union.meryl -readmers readdbmeryl/all_count/union.meryl -peak ~{peak} -vcf vcf.vcf -output out.dump.gz
    >>>
    output {
        File merfin_output = "out.dump.gz"
    }
    runtime {
        docker: "quay.io/chai/merfin:latest"
        cpu: 4
        memory: "32 GB"
        disk: "/tmp 1000 SSD"
    }
    parameter_meta {
    readdb_meryl: {
        description: "meryl intermediate files tar ball for raw reads",
        patterns: ["meryl_files.tar"],
        stream: true,
        localization_optional: true
    }
    seq_meryl: {
        description: "meryl intermediate files tar ball for sequence",
        patterns: ["meryl_files.tar"],
        stream: true,
        localization_optional: true
    }
    sequence_fasta: {
        description: "assembly fasta gzip file",
        patterns: ["*.fasta.gz","*.fa.gz"],
        stream: true,
        localization_optional: true
    }
    vcf: {
        description: "variant VCF file",
        patterns: ["*.vcf.gz"],
        stream: true,
        localization_optional: true
    }
    }
}
