version 1.0

task pretext{
    input {
        File mapbam
    }
    command <<<
        PretextSnapshot -m map.pretext --sequences ~{mapbam} "=full"
    >>>
    output {
        Array[File] pretext_output = glob("*")
    }
    runtime {
        docker: "quay.io/biocontainers/pretext-suite:0.0.1--0"
    }
    parameter_meta {
    mapbam: {
        description: "mapped BAM from Arima",
        extension: ".bam"
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