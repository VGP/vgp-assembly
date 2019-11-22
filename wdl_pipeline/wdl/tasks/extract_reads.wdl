workflow runExtractReads {
    Array[File] inputFiles
    String dockerRepository="tpesout"
    String dockerTag="latest"

    scatter (file in inputFiles) {
        call extractReads {
            input:
                readFile=file,
                dockerRepository=dockerRepository,
                dockerTag=dockerTag
            }
    }

    output {
        Array[File] reads = extractReads.outputFile
    }
}

task extractReads {
    File readFile
    String dockerRepository="tpesout"
    String dockerTag="latest"

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

        export NAME=`basename ${readFile}`
        if [[ $NAME == *bam ]] ; then
            export OUTPUT=`echo $NAME | sed 's/bam$/fastq/'`
            samtools fastq ${readFile} >$OUTPUT
        elif [[ $NAME == *.gz ]] ; then
            export OUTPUT=`echo $NAME | sed 's/.gz$//'`
            gunzip -k -c ${readFile} >$OUTPUT
        elif [[ $NAME == *.fastq ]] || [[ $NAME == *.fq ]] || [[ $NAME == *.fasta ]] || [[ $NAME == *.fa ]] ; then
            export OUTPUT=${readFile}
        else
            echo "Unsupported file type: $NAME"
            exit 1
        fi

        echo $OUTPUT >output
	>>>
	output {
		String outputFilename = read_string("output")
		File outputFile = outputFilename
	}
    runtime {
        cpu: 1
        docker: dockerRepository+"/vgp_base:"+dockerTag
    }
}
