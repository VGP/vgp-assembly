version 1.0

task samtools_fai{
    input {
        File ref
    }
    command <<<
        set -e -x -o pipefail
        ref_name=$(basename ~{ref})
        ref_name=${ref_name%.gz}
        zcat ~{ref} > $ref_name
        samtools faidx $ref_name
        ls
    >>>
    output {
        File fai = glob('*.fai')[0]
    }
    runtime {
        docker: "biocontainers/samtools:v1.9-4-deb_cv1"
        cpu: 4
        memory: "10 GB"
        disk: "/tmp 500 SSD"
    }
    parameter_meta {
    fai: {
        description: "reference index",
        patterns: ["*.fai"]
    }
    }
}
