version 1.0

task merfin_no_header_correction{
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
        mv readdbmeryl/all_count/union_meryl readdbmeryl/all_count/union.meryl
        mv seqdbmeryl/all_count/union_meryl seqdbmeryl/all_count/union.meryl
        zcat ~{sequence_fasta} > sequence.fasta
        zcat ~{vcf} > vcf.vcf
        working_dir=$(pwd)
        #cd /merfin/scripts/reformat_arrow/
        #chmod 777 /merfin/scripts/reformat_arrow/reshape_arrow.sh
        #export merfin="/merfin/scripts/reformat_arrow"
        #/merfin/scripts/reformat_arrow/reshape_arrow.sh "$working_dir"/vcf.vcf
        #ls -lh
        #cd "$working_dir"
        #zcat < vcf.reshaped.vcf.gz > vcf.reshaped.vcf
        ls -lh
        /merfin/build/bin/merfin -vmer -memory1 30 -memory2 150 -comb 15 -sequence sequence.fasta -seqmers seqdbmeryl/all_count/union.meryl -readmers readdbmeryl/all_count/union.meryl -peak ~{peak} -vcf vcf.vcf -output out
        bgzip -c out.polish.vcf > out.polish.vcf.gz
        tabix -p vcf out.polish.vcf.gz
        bcftools consensus out.polish.vcf.gz -f sequence.fasta > out.fasta
        gzip out.fasta
    >>>
    output {
        File merfin_vcf = "out.polish.vcf"
        File merfin_debug = "out.debug"
        File merfin_genome = "out.fasta.gz"

    }
    runtime {
        docker: "quay.io/chai/merfin:latest"
        cpu: 8
        memory: "200 GB"
        disk: "/tmp 620 SSD"
    }
    parameter_meta {
    readdb_meryl: {
        description: "meryl intermediate files tar ball for raw reads",
        patterns: ["meryl_files.tar"]
    }
    seq_meryl: {
        description: "meryl intermediate files tar ball for sequence",
        patterns: ["meryl_files.tar"]
    }
    sequence_fasta: {
        description: "assembly fasta gzip file",
        patterns: ["*.fasta.gz","*.fa.gz"]
    }
    vcf: {
        description: "variant VCF file",
        patterns: ["*.vcf.gz"]
    }
    }
}
