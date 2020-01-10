version 1.0

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
