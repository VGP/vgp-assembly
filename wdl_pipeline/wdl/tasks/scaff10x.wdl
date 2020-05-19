version 1.0

workflow runScaff10x {
	call scaff10x
}

task scaff10x {
    input {
        File assemblyFasta
        Array[File] readFiles10x
        String sampleName
        Int threadCount
        Int? memoryGigabyte=32
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

        # link assembly
        ln -s ~{assemblyFasta} asm.fasta

        # link input files
        mkdir input
        for RF in ~{sep=" " readFiles10x} ; do
            if [[ ! $RF == *.fastq.gz ]] ; then
                echo "10x read files expected to be in .fastq.gz format: $RF"
                exit 1
            fi
            ln -s $RF input/$(basename $RF)
        done

        # make input.dat
        for RF in $(ls input/) ; do
            # sort into R1 and R2
            if [[ $RF == *R1* ]] && [[ $RF == *R2* ]] ; then
                echo "Cannot determine R1 or R2 from filename (found both): $RF"
                exit 1
            elif [[ $RF == *R1* ]] ; then
                echo "q1=$(pwd)/input/$RF" >>input.dat
            elif [[ $RF == *R2* ]] ; then
                echo "q2=$(pwd)/input/$RF" >>input.dat
            else
                echo "Cannot determine R1 or R2 from filename (found neither): $RF"
                exit 1
            fi
        done

        # validate input.dat
        if [[ ! -s input.dat ]] ; then
            echo "Check input.dat"
            exit -1
        fi

        # run script
        export tools=/root/tools
        export SLURM_CPUS_PER_TASK=~{threadCount}
        bash /root/scripts/scaff10x/scaff10x_v4.1.sh ~{sampleName}
	>>>
	output {
		File scaffoldedFasta = sampleName + ".scaff10x.fasta"
	}
    runtime {
        cpu: 32
        memory: "72 GB"
        docker: dockerImage
    }
}
