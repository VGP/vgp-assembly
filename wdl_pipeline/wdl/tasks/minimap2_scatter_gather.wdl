version 1.0

workflow runMinimap2ScatterGather {
  input {
    File refFasta
    Array[File] readFiles
    Int threadCount
    String? minimapPreset
    String? samtoolsFilter
    String dockerRepository
    String dockerTag
  }

  String dockerImage = "${dockerRepository}/vgp_minimap2:${dockerTag}"

	call minimap2_idx as idx {
    input:
      refFasta=refFasta,
      threadCount=threadCount,
      minimapPreset=minimapPreset,
      dockerImage=dockerImage
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
        dockerImage=dockerImage
      }
    }

	call minimap2_merge as mrg {
    input:
      outputBase=idx.outputBase,
      alignFiles=aln.minimap2Bam,
      alignIndexFiles=aln.minimap2BamIdx,
      threadCount=threadCount,
      dockerImage=dockerImage
	}

  output {
    File alignment = mrg.minimap2Bam
    File alignmentIdx = mrg.minimap2BamIdx
  }
}

task minimap2_idx {
  input {
    File refFasta
    Int threadCount
    String? minimapPreset
    String dockerImage
  }

  String defaultMinimapPreset = select_first([minimapPreset, ""])

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
    filename=$(basename -- "~{refFasta}")
    prefix="${filename%.*}"
    suffix="${filename##*.}"
    if [[ "${suffix}" == "gz" ]]; then
      prefix="${prefix##*.}"
    fi

    # index ref
    export SLURM_CPUS_PER_TASK=~{threadCount}
    bash /root/scripts/minimap2/minimap2_idx.sh ~{refFasta} ~{defaultMinimapPreset}

    mkdir output
    mv ${prefix}.idx output
	>>>

	output {
	    String outputBase = read_string("outputBase")
		File minimap2Index = glob("output/*.idx")[0]
	}

  runtime {
    docker: dockerImage
    cpu: threadCount
    memory: "32 GB"
  }

  parameter_meta {
    refFasta: {
      description: "Reference genome in FASTA format",
      stream: true,
      localization_optional: true
    }
  }
}

task minimap2_align {
  input {
    File refFasta
    File refFastaIdx
    File readFile
    Int threadCount
    String? minimapPreset
    String? samtoolsFilter
    String dockerImage
  }

  String defaultMinimapPreset = select_first([minimapPreset, ""])
  String defaultSamtoolsFilter = select_first([samtoolsFilter, ""])

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
    filename=$(basename -- "~{readFile}")
    prefix="${filename%.*}"
    suffix="${filename##*.}"
    if [[ "${suffix}" == "gz" ]]; then
      prefix="${prefix##*.}"
    fi

    ln -s ~{refFasta}
    ln -s ~{refFastaIdx}
    ln -s ~{readFile}
    echo `basename ~{readFile}` > input.fofn

    export SLURM_CPUS_PER_TASK=~{threadCount}
    bash /root/scripts/minimap2/minimap2.sh \
      ${filename} 1 "~{minimapPreset}" "~{samtoolsFilter}"

    mkdir output
    mv *.bam *.bai output
	>>>

	output {
		File minimap2Bam = glob("output/*.bam")[0]
		File minimap2BamIdx = glob("output/*.bai")[0]
	}

  runtime {
    docker: dockerImage
    cpu: threadCount
    memory: "32 GB"
  }

  parameter_meta {
    refFasta: {
      description: "Reference FASTA file",
      stream: true,
      localization_optional: true
    }
    refFastaIdx: {
      description: "Reference Minimap2 Index file",
      stream: true,
      localization_optional: true
    }
    readFile: {
      description: "Read file",
      stream: true,
      localization_optional: true
    }
  }
}

task minimap2_merge {
  input {
    String outputBase
    Array[File] alignFiles
    Array[File] alignIndexFiles
    Int threadCount
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

    for RF in ~{sep=" " alignFiles} ; do
      ln -s $RF
    done

    for RF in ~{sep=" " alignIndexFiles} ; do
      ln -s $RF
    done

    export SLURM_CPUS_PER_TASK=~{threadCount}
    bash /root/scripts/minimap2/merge.sh ~{outputBase}
	>>>

	output {
		File minimap2Bam = "${outputBase}.bam"
		File minimap2BamIdx = "${outputBase}.bam.bai"
	}

  runtime {
    docker: dockerImage
    cpu: threadCount
    memory: "32 GB"
  }

  parameter_meta {
    alignFiles: {
      description: "BAM files to merge",
      stream: true,
      localization_optional: true
    }
    alignIndexFiles: {
      description: "BAM index files",
      stream: true,
      localization_optional: true
    }
  }
}
