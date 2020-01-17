version 1.0


workflow PolishAssembly {
    input {
        File ASSEMBLY_FILE
        Array[File] READ_FILES
        String SAMPLE_NAME
        File MARGIN_POLISH_PARAMS
        String? MINIMAP_PRESET
        String? SAMTOOLS_FILTER
        Int THREAD_COUNT
        Int? MARGINPOLISH_MEMORY_GB
        String? DOCKER_REPOSITORY
        String? DOCKER_TAG
    }

    String defaultMinimapPreset = select_first([MINIMAP_PRESET, "map-ont"])
    String defaultSamtoolsFilter = select_first([SAMTOOLS_FILTER, "-F 0x904"])
    Int defaultMarginPolishMemoryGB = select_first([MARGINPOLISH_MEMORY_GB, 8])
    String defaultDockerRepository = select_first([DOCKER_REPOSITORY, "tpesout"])
    String defaultDockerTag = select_first([DOCKER_TAG, "latest"])

    String mm2DockerImage = "${defaultDockerRepository}/vgp_minimap2:${defaultDockerTag}"

  scatter (file in READ_FILES) {
    call extractReads as xtract {
      input:
        readFile=file,
        dockerRepository=defaultDockerRepository,
        dockerTag=defaultDockerTag
      }
  }

	call minimap2_idx as idx {
    input:
      refFasta=ASSEMBLY_FILE,
      threadCount=THREAD_COUNT,
      minimapPreset=defaultMinimapPreset,
      dockerImage=mm2DockerImage
	}

  scatter (file in xtract.outputFile) {
    call minimap2_align as aln {
      input:
        refFasta=ASSEMBLY_FILE,
        refFastaIdx=idx.minimap2Index,
        readFile=file,
        threadCount=THREAD_COUNT,
        minimapPreset=defaultMinimapPreset,
        samtoolsFilter=defaultSamtoolsFilter,
        dockerImage=mm2DockerImage
      }
    }

	call minimap2_merge as mrg {
    input:
      outputBase=idx.outputBase,
      alignFiles=aln.minimap2Bam,
      alignIndexFiles=aln.minimap2BamIdx,
      threadCount=THREAD_COUNT,
      dockerImage=mm2DockerImage
	}

    call marginPolish as polish {
        input:
            sampleName=SAMPLE_NAME,
            alignmentBam=mrg.minimap2Bam,
            alignmentBamIdx=mrg.minimap2BamIdx,
            referenceFasta=ASSEMBLY_FILE,
            parameters=MARGIN_POLISH_PARAMS,
            featureType="",
            threadCount=THREAD_COUNT,
            memoryGigabyte=defaultMarginPolishMemoryGB,
            dockerRepository=defaultDockerRepository,
            dockerTag=defaultDockerTag
    }

    output {
        File polishedAssembly = polish.polishedFasta
    }

}



task extractReads {
  input {
    File readFile
    String dockerRepository
    String dockerTag
  }

  String dockerImage = "${dockerRepository}/vgp_base:${dockerTag}"

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

    module load samtools

    filename=$(basename -- "~{readFile}")
    prefix="${filename%.*}"
    suffix="${filename##*.}"

    mkdir output

    if [[ "$suffix" == "bam" ]] ; then
      samtools fastq ~{readFile} > output/${prefix}.fq
    elif [[ "$suffix" == "gz" ]] ; then
      gunzip -k -c ~{readFile} > output/${prefix}
    elif [[ "$suffix" != "fastq" ]] && [[ "$suffix" != "fq" ]] && [[ "$suffix" != "fasta" ]] && [[ "$suffix" != "fa" ]] ; then
      echo "Unsupported file type: ${suffix}"
      exit 1
    fi
	>>>

	output {
	  File outputFile = flatten([glob("output/*"), [readFile]])[0]
	}

  runtime {
    docker: dockerImage
    cpu: 1
  }

  parameter_meta {
    readFile: {
      description: "Reads file in BAM, FASTQ, or FASTA format (optionally gzipped)",
      stream: true,
      localization_optional: true
    }
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
    ln -s ~{refFasta}
    filename=$(basename -- "~{refFasta}")
    prefix="${filename%.*}"
    suffix="${filename##*.}"
    if [[ "${suffix}" == "gz" ]]; then
      prefix="${prefix##*.}"
    fi
    echo $prefix >outputBase

    # index ref
    export SLURM_CPUS_PER_TASK=~{threadCount}
    bash /root/scripts/minimap2/minimap2_idx.sh $filename ~{defaultMinimapPreset}

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
    filename=$(basename -- "~{refFasta}")
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

task marginPolish {
  input {
    String sampleName
    File alignmentBam
    File alignmentBamIdx
    File referenceFasta
    File parameters
    String? featureType
    Int threadCount
    Int memoryGigabyte
    String dockerRepository
    String dockerTag
  }

  String defaultFeatureType = select_first([featureType, ""])
  String dockerImage = "${dockerRepository}/vgp_marginpolish:${dockerTag}"

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


    export FEATURE_PARAM=""
    if [[ ! -z "~{defaultFeatureType}" ]] ; then
        export FEATURE_PARAM="-F ~{defaultFeatureType}"
    fi

    mkdir output
    ln -s ~{alignmentBam}
    ln -s ~{alignmentBamIdx}
    module load marginPolish

    marginPolish \
      `basename ~{alignmentBam}` \
      ~{referenceFasta} \
      ~{parameters} \
      -t ~{threadCount} \
      -o output/~{sampleName}.marginPolish \
      $FEATURE_PARAM
	>>>

	output {
		File polishedFasta = "output/${sampleName}.marginPolish.fa"
		Array[File] helenImages = glob("output/*.h5")
	}

  runtime {
    docker: dockerImage
    cpu: threadCount
    memory: "${memoryGigabyte} GB"
  }

  parameter_meta {
    alignmentBam: {
      description: "BAM file to polish",
      stream: true,
      localization_optional: true
    }
    alignmentBamIdx: {
      description: "BAM index file",
      stream: true,
      localization_optional: true
    }
    referenceFasta: {
      description: "Reference genome in FASTA format",
      stream: true,
      localization_optional: true
    }
  }
}

