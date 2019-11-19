workflow runStats {
	call stats
}

task stats {
    File assemblyFasta
    Int expectedGenomeSize
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

        export VGP_PIPELINE=/root/scripts

        ln -s ${assemblyFasta}
        export REF=`basename ${assemblyFasta}`
        echo $REF | sed 's/.fasta//' | sed 's/.fa//' >outputBase

        bash /root/scripts/stats/asm_stats.sh $REF ${expectedGenomeSize}

        tar czvf $(cat outputBase).stats.tar.gz *.stats
	>>>
	output {
		String outputBase = read_string("outputBase")
		File statsTarball = outputBase + ".stats.tar.gz"
	}
    runtime {
        cpu: 1
        memory: 1 + " GB"
        docker: dockerRepository+"/vgp_stats:"+dockerTag
    }
}
