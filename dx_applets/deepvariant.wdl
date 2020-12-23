version 1.0

task deepvariant{
    input {
        File bam
        File bai
        File ref
        File fai
        String model='WGS'
    }
    command <<<
        set -e -x -o pipefail
        cat ~{bam} > input.bam
        cat ~{bai} > input.bam.bai
        read_name=$(basename ~{bam} .bam)
        ref_name=$(basename ~{ref})
        ref_name=${ref_name%.gz}
        fai_name=$(basename ~{fai})
        zcat ~{ref} > $ref_name
        cat ~{fai} > $fai_name
        /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$ref_name --reads=input.bam --output_vcf=$read_name.vcf.gz --num_shards=$(nproc)
    >>>
    output {
        File vcf = glob('*.vcf.gz')[0]
        File tbi = glob('*.vcf.gz.tbi')[0]
        File html = glob('*.html')[0]
    }
    runtime {
        docker: "google/deepvariant:1.1.0"
        cpu: 26
        memory: "24 GB"
        disk: "/tmp 1000 SSD"
    }
    parameter_meta {
    bam: {
        description: "input bam file",
        patterns: ["*.bam"]
    }
    bai: {
        description: "input bai file",
        patterns: ["*.bai"]
    }
    ref: {
        description: "reference",
        patterns: ["*.fasta.gz","*.fa.gz"]
    }
    fai: {
        description: "reference index",
        patterns: ["*.fai"]
    }
    }
}
