workflow helloScaff10x {
	call scaff10x
}

task scaff10x {
    File assemblyFasta
    Array[File] readFiles10x
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

        # get name of output file
        REF=`basename ${assemblyFasta} | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

        # link input files
        mkdir input
        for RF in ${sep=" " readFiles10x} ; do
            ln -s $RF input/$(basename $RF)
        done

        # make and validate input.dat
        ls input/*_R[1-2]_001.fastq.gz | awk '{if (NR%2==1) {print "q1="$0} else {print "q2="$0}}' > input.dat
        if [[ ! -s input.dat ]] ; then
            echo "Check input.dat"
            exit -1
        fi

        # run script
        export tools=/root/tools
        bash /root/scripts/scaff10x/scaff10x_v4.1.sh

        mv $REF.scaff10x.fasta ${sampleName}.scaff10x.fasta
	>>>
	output {
		File scaffoldedFasta = sampleName + ".scaff10x.fasta"
	}
    runtime {
        cpu: 32
        memory: "72 GB"
        docker: "tpesout/vgp_scaff10x"
    }
}
