workflow runPurgeDups {
	call purge_dups
}

task purge_dups {
    File assemblyFasta
    Array[File] readFiles
    String minimapPreset
    String sampleName
    Int threadCount
    Int memoryGigabyte
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

        # link ref and get name of output file
        ASM=`basename ${assemblyFasta}`
        ln -s ${assemblyFasta}

        # link read files and generate input.fofn
        for RF in ${sep=" " readFiles} ; do
            ln -s $RF
            echo `basename $RF` >>input.fofn
        done

        # environment variables for scripts
        export SLURM_CPUS_PER_TASK=${threadCount}
        export tools="/root/tools"

        # index ref (if not present)
        bash /root/scripts/purge_dups/minimap2_idx.sh $ASM ${minimapPreset}

        # get read PAFs
        IDX=1
        for RF in ${sep=" " readFiles} ; do
            bash /root/scripts/purge_dups/minimap2.sh $ASM $IDX ${minimapPreset}
            IDX=$(($IDX + 1))
        done

        # get self align
        bash /root/scripts/purge_dups/minimap2_self.sh $ASM

        # run purge dups
        bash /root/scripts/purge_dups/purge_dups.sh $ASM

        # cleanup
        mv purged.fa ${sampleName}.purged.fa
        mv hap.fa ${sampleName}.hap.fa

	>>>
	output {
		File primary = sampleName + ".purged.fa"
		File alternate = sampleName + ".hap.fa"
	}
    runtime {
        cpu: threadCount
        memory: memoryGigabyte + " GB"
        docker: dockerRepository+"/vgp_purge_dups:"+dockerTag
    }
}
