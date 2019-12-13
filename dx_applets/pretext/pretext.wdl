version 1.0

task pretext{
    input {
        File mapbam
    }
#    Int disk_space = ceil(2.0 * sum(size(mapbam, "G")))
    command <<<
        set -e -x -o pipefail
        cat ~{mapbam} | samtools view -h - | PretextMap -o map.pretext
        PretextSnapshot -m map.pretext --sequences "=full"
        ls
    >>>
    output {
        Array[File] pretext_output = glob("map_snapshots/map_FullMap.png")
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
        Array[File] pretext_output = pretext.pretext_output
    }
}