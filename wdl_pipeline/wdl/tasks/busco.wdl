workflow runBusco {
	call busco
}

task busco {
    File assemblyFasta
    Int threadCount
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

        # get name of output file
        ASSEMBLY=`basename ${assemblyFasta} | sed 's/.fasta$//g' | sed 's/.fa$//g' `
        echo $ASSEMBLY >outputBase

        ln -s ${assemblyFasta}

        export SLURM_CPUS_PER_TASK=${threadCount}
        export tools=/root/tools

        bash /root/scripts/busco/busco.sh `basename ${assemblyFasta}`

        tar czvf $ASSEMBLY.busco.tar.gz run_$ASSEMBLY

	>>>
	output {
		String outputBase = read_string("outputBase")
		File outputTarball = outputBase + ".busco.tar.gz"
	}
    runtime {
        cpu: threadCount
        memory: "42 GB"
        docker: dockerRepository+"/vgp_busco:"+dockerTag
    }
}
