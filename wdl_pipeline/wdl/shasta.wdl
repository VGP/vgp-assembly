workflow helloShasta {
	call shasta
}

task shasta {
    Array[File] readFiles
    String sampleName
    Int threadCount
    Int memoryGigabyte
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

        module load shasta
        shasta --input ${sep=" --input " readFiles} --threads ${threadCount}
        mv ShastaRun/Assembly.fasta ${sampleName}.shasta.fasta
        mv ShastaRun ${sampleName}
        tar czvf ${sampleName}.shasta.tar.gz ${sampleName}/*
	>>>
	output {
		File assemblyFasta = sampleName + ".shasta.fasta"
		File assemblyExtra = sampleName + ".shasta.tar.gz"
	}
    runtime {
        cpu: threadCount
        memory: memoryGigabyte + " GB"
        docker: "tpesout/vgp_shasta"
    }
}
