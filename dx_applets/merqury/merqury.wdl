version 1.0

task merqury{
    input {
        File readdb_meryl
        File asm1_fasta
        File? mat_meryl
        File? pat_meryl
        File? asm2_fasta
        String output_prefix
    }
    command <<<
        set -e -x -o pipefail
        cmd=(merqury.sh ~{readdb_meryl})
        if [[ -f "~{mat_meryl}" ]]; then
            cmd+=(~{mat_meryl})
        fi
        if [[ -f "~{pat_meryl}" ]]; then
            cmd+=(~{pat_meryl})
        fi
        cmd+=(~{asm1_fasta})
        if [[ -f "~{asm2_fasta}" ]]; then
            cmd+=(~{asm2_fasta})
        fi
        cmd+=~{output_prefix}
        ls
    >>>
    output {
        File pretext_snapshot = "map_snapshots/map_FullMap.png"
        File map_pretext = "map.pretext"
    }
    runtime {
        docker: "quay.io/chai/pretext:0.0.2"
        cpu: 4
        memory: "50 GB"
#        disk: "/tmp ${disk_space} SSD"
        disk: "/tmp 50 SSD"
    }
    parameter_meta {
    mapbam: {
        description: "mapped BAM from Arima",
        extension: ".bam",
        stream: true,
        localization_optional: true
    }
    }
}

workflow pretext_workflow{
    input {
        File mapbam
    }
    call pretext {
        input:
            mapbam = mapbam
    }
    output {
        File pretext_snapshot = pretext.pretext_snapshot
        File map_pretext = pretext.map_pretext
    }
}