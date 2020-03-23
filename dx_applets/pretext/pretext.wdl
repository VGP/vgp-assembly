version 1.0

task pretext{
    input {
        File mapbam
        Int? color
    }
    Int actual_color = select_first([color,5])
#    Int disk_space = ceil(2.0 * sum(size(mapbam, "G")))
    command <<<
        set -e -x -o pipefail
        cat ~{mapbam} | samtools view -h - | PretextMap -o map.pretext
        PretextSnapshot --colourMap ~{actual_color} -m map.pretext --sequences "=full"
        ls
    >>>
    output {
        File pretext_snapshot = "map_snapshots/map_FullMap.png"
        File map_pretext = "map.pretext"
    }
    runtime {
        docker: "quay.io/chai/pretext:0.0.1"
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