version 1.0

task qv_assessment{
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
        cmd=(bash /merqury/eval/qv.sh all_count/union.meryl)
        cmd+=(assem1.fa)
        if [[ -f "~{asm2_fasta}" ]]; then
            zcat ~{asm2_fasta} >  assem2.fa
            cmd+=(assem2.fa)
        fi
        cmd+=(~{output_prefix})
        "${cmd[@]}"
    >>>
    output {
        File qv = "${output_prefix}.qv"
    }
    runtime {
        docker: "quay.io/chai/merqury:0.0.1"
        cpu: 4
        memory: "32 GB"
        disk: "/tmp 1000 SSD"
    }
    parameter_meta {
    readdb_meryl: {
        description: "meryl intermediate files tar ball",
        extension: ".tar",
        stream: true,
        localization_optional: true
    }
    asm1_fasta: {
        description: "first assembly",
        extension: ".gz",
        stream: true,
        localization_optional: true
    }
    asm2_fasta: {
        description: "second assembly",
        extension: ".gz",
        stream: true,
        localization_optional: true
    }
    }
}

workflow qv_workflow{
    input {
        File readdb_meryl
        File asm1_fasta
        File? asm2_fasta
    }
    call qv_assessment {
        input:
            readdb_meryl = readdb_meryl,
            asm1_fasta = asm1_fasta,
            asm2_fasta = asm2_fasta
    }
    output {
        File qv = qv_assessment.qv
    }
}