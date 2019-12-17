workflow runMarginPolish {
	call marginPolish
}

task marginPolish {
    String sampleName
    File alignmentBam
    File alignmentBamIdx
    File referenceFasta
    File parameters
    String featureType=""
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

        export FEATURE_PARAM="-F ${featureType}"
        if [ -z "${featureType}" ]; then
            export FEATURE_PARAM=""
        fi
        mkdir output
        ln -s ${alignmentBam}
        ln -s ${alignmentBamIdx}
        module load marginPolish

        marginPolish \
            `basename ${alignmentBam}` \
            ${referenceFasta} \
            ${parameters} \
            -t ${threadCount} \
            -o output/${sampleName}.marginPolish \
            $FEATURE_PARAM
	>>>
	output {
		File polishedFasta = "output/" + sampleName + ".marginPolish.fa"
		Array[File] helenImages = glob("output/*.h5")
	}
    runtime {
        cpu: threadCount
        memory: memoryGigabyte + " GB"
        docker: dockerRepository+"/vgp_marginpolish:"+dockerTag
    }
}
