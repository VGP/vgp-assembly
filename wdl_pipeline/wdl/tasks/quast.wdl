version 1.0

workflow runQuast {
	call quast
}

task quast {
    input {
        File assemblyFasta
        File referenceFasta
        String sampleName
        Int threadCount
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

        # link input files
        ln -s ~{assemblyFasta}
        ln -s ~{referenceFasta}

        module load python/2.7
        python /root/tools/quast/quast-5.0.2/quast-lg.py \
            -t ~{threadCount} \
            -r `basename ~{referenceFasta}` \
            --min-identity 80 \
            --fragmented \
            --large \
            -o ~{sampleName}.quast \
            `basename ~{assemblyFasta}`

        tar czvf ~{sampleName}.quast.tar.gz ~{sampleName}.quast

	>>>
	output {
		File outputTarball = sampleName + ".quast.tar.gz"
	}
    runtime {
        cpu: threadCount
        memory: "8 GB"
        docker: dockerImage
    }
}
