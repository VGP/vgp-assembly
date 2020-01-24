version 1.0

workflow runExtractReads {
    input {
        Array[File] inputFiles
        String dockerImage
    }

    scatter (file in inputFiles) {
        call extractReads {
            input:
                readFile=file,
                dockerImage=dockerImage
        }
    }

    output {
        Array[File] reads = extractReads.outputFile
    }
}

task extractReads {
    input {
        File readFile
        String dockerImage
    }


    command <<<
        # initialize modules
        source /usr/local/Modules/init/bash
        module use /root/modules/
        # Set the exit code of a pipeline to that of the rightmost command
        # to exit with a non-zero status, or zero if all commands of the pipeline exit
        set -o pipefail
        # cause a bash script to exit immediately when a command fails
        set -e
        # cause the bash shell to treat unset variables as an error and exit immediately
        set -u
        # echo each line of the script to stdout so we can see what is happening
        # to turn off echo do 'set +o xtrace'
        set -o xtrace

        module load samtools

        filename=$(basename -- "~{readFile}")
        prefix="${filename%.*}"
        suffix="${filename##*.}"

        mkdir output

        if [[ "$suffix" == "bam" ]] ; then
            samtools fastq ~{readFile} > output/${prefix}.fq
        elif [[ "$suffix" == "gz" ]] ; then
            gunzip -k -c ~{readFile} > output/${prefix}
        elif [[ "$suffix" != "fastq" ]] && [[ "$suffix" != "fq" ]] && [[ "$suffix" != "fasta" ]] && [[ "$suffix" != "fa" ]] ; then
            echo "Unsupported file type: ${suffix}"
            exit 1
        fi
    >>>

    output {
        File outputFile = flatten([glob("output/*"), [readFile]])[0]
    }

    runtime {
#    String dockerImage = "${dockerRepository}/vgp_base:${dockerTag}"
        docker: dockerImage
        cpu: 1
    }

    parameter_meta {
        readFile: {description: "Reads file in BAM, FASTQ, or FASTA format (optionally gzipped)"}
    }
}
