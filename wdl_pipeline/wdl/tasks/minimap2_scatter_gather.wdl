workflow runMinimap2ScatterGather {

    File refFasta
    Array[File] readFiles
    Int threadCount
    String minimapPreset=""
    String samtoolsFilter=""
    String dockerRepository="tpesout"
    String dockerTag="latest"

	call minimap2_idx as idx {
        input:
            refFasta=refFasta,
            threadCount=threadCount,
            minimapPreset=minimapPreset,
            dockerRepository=dockerRepository,
            dockerTag=dockerTag,
	}

    scatter (file in readFiles) {
        call minimap2_align as aln {
            input:
                refFasta=refFasta,
                refFastaIdx=idx.minimap2Index,
                readFile=file,
                threadCount=threadCount,
                minimapPreset=minimapPreset,
                samtoolsFilter=samtoolsFilter,
                dockerRepository=dockerRepository,
                dockerTag=dockerTag
            }
    }

	call minimap2_merge as mrg {
        input:
            outputBase=idx.outputBase,
            alignFiles=aln.minimap2Bam,
            alignIndexFiles=aln.minimap2BamIdx,
            threadCount=threadCount,
            dockerRepository=dockerRepository,
            dockerTag=dockerTag,
	}

    output {
        File alignment = mrg.minimap2Bam
        File alignmentIdx = mrg.minimap2BamIdx
    }

}

task minimap2_idx {
    File refFasta
    Int threadCount
    String minimapPreset=""
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

        # index ref
        export SLURM_CPUS_PER_TASK=${threadCount}
        bash /root/scripts/minimap2/minimap2_idx.sh $REF ${minimapPreset}

	>>>
	output {
		String outputBase = read_string("outputBase")
		File minimap2Index = outputBase + ".idx"
	}
    runtime {
        cpu: threadCount
        memory: "32 GB"
        docker: dockerRepository+"/vgp_minimap2:"+dockerTag
    }
}
task minimap2_align {
    File refFasta
    File refFastaIdx
    Int threadCount
    File readFile
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
        READ_ID=`basename ${readFile}`
        READ_ID=`echo $READ_ID | sed 's/.fasta.$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
        READ_ID=`echo $READ_ID | sed 's/.fastq.$//g' | sed 's/.fastq$//g' | sed 's/.fq$//g' | sed 's/.fastq.gz$//g' | sed 's/.fq.gz$//g'`

        echo $READ_ID > outputBase

        ln -s ${refFasta}
        ln -s ${refFastaIdx}
        ln -s ${readFile}
        echo `basename ${readFile}` >input.fofn

        export SLURM_CPUS_PER_TASK=${threadCount}

        bash /root/scripts/minimap2/minimap2.sh `basename ${refFasta}` 1 "${minimapPreset}" "${samtoolsFilter}"

	>>>
	output {
		String outputBase = read_string("outputBase")
		File minimap2Bam = outputBase + ".sort.bam"
		File minimap2BamIdx = outputBase + ".sort.bam.bai"
	}
    runtime {
        cpu: threadCount
        memory: "32 GB"
        docker: dockerRepository+"/vgp_minimap2:"+dockerTag
    }
}


task minimap2_merge {
    String outputBase
    Array[File] alignFiles
    Array[File] alignIndexFiles
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

        for RF in ${sep=" " alignFiles} ; do
            ln -s $RF
        done

        for RF in ${sep=" " alignIndexFiles} ; do
            ln -s $RF
        done

        export SLURM_CPUS_PER_TASK=${threadCount}

        bash /root/scripts/minimap2/merge.sh ${outputBase}

	>>>
	output {
		File minimap2Bam = outputBase + ".bam"
		File minimap2BamIdx = outputBase + ".bam.bai"
	}
    runtime {
        cpu: threadCount
        memory: "32 GB"
        docker: dockerRepository+"/vgp_minimap2:"+dockerTag
    }
}