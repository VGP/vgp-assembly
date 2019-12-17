workflow runMinimap2 {
	call minimap2 
}

task minimap2 {
    File refFasta
    Array[File] readFiles
    Int threadCount
    String minimapPreset=""
    String samtoolsFilter=""
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
        REF=`basename ${refFasta}`
        echo $REF | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' >outputBase

        ln -s ${refFasta}

        for RF in ${sep=" " readFiles} ; do
            ln -s $RF
            echo `basename $RF` >>input.fofn
        done

        export SLURM_CPUS_PER_TASK=16

        # index ref (if not present)
        bash /root/scripts/minimap2/minimap2_idx.sh $REF ${minimapPreset}

        # align
        IDX=1
        for RF in ${sep=" " readFiles} ; do
            bash /root/scripts/minimap2/minimap2.sh $REF $IDX ${minimapPreset} "${samtoolsFilter}"
            IDX=$(($IDX + 1))
        done

        bash /root/scripts/minimap2/merge.sh `cat outputBase`

	>>>
	output {
		String outputBase = read_string("outputBase")
		File minimap2Bam = outputBase + ".bam"
		File minimap2BamIdx = outputBase + ".bam.bai"
	}
    runtime {
        cpu: threadCount
        memory: "32 GB"
        docker: dockerRepository+"/vgp_minimap2:"+dockerTag
    }
}
