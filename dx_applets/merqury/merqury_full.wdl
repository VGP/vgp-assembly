version 1.0

task merqury{
    input {
        File readdb_meryl
        File asm1_fasta
        File? asm2_fasta
        String output_prefix
    }
    command <<<
        set -e -x -o pipefail
        echo MERQURY
        echo $MERQURY
        ls $MERQURY
        tar -xvf ~{readdb_meryl}
        zcat ~{asm1_fasta} > assem1.fa
        mv all_count/union_meryl all_count/union.meryl
        cmd=(bash /merqury/merqury.sh all_count/union.meryl)
        cmd+=(assem1.fa)
        if [[ -f "~{asm2_fasta}" ]]; then
            zcat ~{asm2_fasta} >  assem2.fa
            cmd+=(assem2.fa)
        fi
        cmd+=(~{output_prefix})
        "${cmd[@]}"
        ls
    >>>
    output {
        Array[File] qv = glob('*.qv')
        Array[File] hist = glob('*hist*')
        Array[File] png = glob('*.png')
        Array[File] stats = glob('*.stats')
        Array[File] log = glob('*.log')
        Array[File] bed = glob('*.bed')
        Array[File] tdf = glob('*.tdf')
        Array[File] filt = glob('*.filt')
    }
    runtime {
        docker: "quay.io/chai/merqury:1.1"
        cpu: 4
        memory: "32 GB"
        disk: "/tmp 1000 SSD"
    }
    parameter_meta {
    readdb_meryl: {
        description: "meryl intermediate files tar ball",
        patterns: ["meryl_files.tar"],
    }
    asm1_fasta: {
        description: "first assembly",
        patterns: ["*.fasta.gz","*.fa.gz"],
    }
    asm2_fasta: {
        description: "second assembly",
        patterns: ["*.fasta.gz","*.fa.gz"],
    }
    }
}
