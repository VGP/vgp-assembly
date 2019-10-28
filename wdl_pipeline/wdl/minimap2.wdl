workflow helloMinimap2 {
	call minimap2 
}

task minimap2 {
    File refFasta
    File readFile
    String dockerMinimap2Tag
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
        REF=`basename ${refFasta}`
        echo $REF | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' >outputBase
        READS=`basename ${readFile}`

        ln -s ${refFasta}
        ln -s ${readFile}
        echo $READS >input.fofn
        export SLURM_CPUS_PER_TASK=16

        # index ref (if not present)
        if [ ! -e $REF.idx ]; then
            bash /root/scripts/minimap2/minimap2_idx.sh $REF
        fi

        # align
        if [ ! -e $REF.bam ]; then
            bash /root/scripts/minimap2/minimap2.sh $REF 1
        fi

        bash /root/scripts/minimap2/merge.sh `cat outputBase`

	>>>
	output {
		String outputBase = read_string("outputBase")
		File minimap2Bam = outputBase + ".bam"
		File minimap2BamIdx = outputBase + ".bam.bai"
	}
    runtime {
        cpu: 16
        memory: "32 GB"
        docker: "tpesout/vgp_minimap2"
    }
}
