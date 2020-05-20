version 1.0

workflow runMerqury {
	call merqury
}

task merqury {
    input {
        File assemblyFasta
        File? altHapFasta
        File kmerTarball
        File? matKmerTarball
        File? patKmerTarball
        String sampleName
        Int threadCount
        String dockerImage
    }

	command <<<
        # don't initialize modules (image was built with different infrastructure)
        #source /usr/local/Modules/init/bash
        #module use /root/modules/
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

        # initilize command
        cmd=(/merqury/merqury.sh)

        # get kmers
        tar xvf ~{kmerTarball}
        cmd+=($(basename ~{kmerTarball} | sed 's/.gz$//' | sed 's/.tar$//'))
        if [[ -f "~{matKmerTarball}" && -f "~{patKmerTarball}" ]]; then
            tar xvf ~{matKmerTarball}
            tar xvf ~{patKmerTarball}
            cmd+=($(basename ~{matKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//'))
            cmd+=($(basename ~{patKmerTarball} | sed 's/.gz$//' | sed 's/.tar$//'))
        fi

        # link input files
        ln -s ~{assemblyFasta} asm.fasta
        cmd+=(asm.fasta)
        if [[ -f "~{altHapFasta}" ]]; then
            ln ~{altHapFasta} altHap.fasta
            cmd+=(altHap.fasta)
        fi

        # prep output
        cmd+=(~{sampleName}.merqury)

        # run command
        echo "${cmd[@]}"
        ${cmd[@]}

        # get output
        tar czvf ~{sampleName}.merqury.tar.gz ~{sampleName}.merqury*

	>>>
	output {
		File outputTarball = sampleName + ".merqury.tar.gz"
	}
    runtime {
        cpu: threadCount
        memory: "8 GB"
        docker: dockerImage
    }
}
