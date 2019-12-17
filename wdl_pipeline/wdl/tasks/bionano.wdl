version 1.0

task bionano_solve {
    input {
        Array[File] bionanoFiles
        File assemblyFasta
        String sampleName
        Int threadCount
        String dockerRepository="tpesout"
        String dockerTag="latest"
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

        # potential files and cmaps
        export BSPQI_CMAP=""
        export BSPQI_BNX=""
        export BSSSI_CMAP=""
        export BSSSI_BNX=""
        export DLE1_CMAP=""
        export DLE1_BNX=""

        # look through all files
        shopt -s nocasematch
        for FILE in ~{sep=" " bionanoFiles} ; do
            if [[ $FILE =~ "BspQI" ]] ; then
                if [[ $FILE =~ "cmap" ]] ; then export BSPQI_CMAP=$FILE ; fi
                if [[ $FILE =~ "bnx" ]] ; then export BSPQI_BNX=$FILE ; fi
            fi
            if [[ $FILE =~ "BssSI" ]] ; then
                if [[ $FILE =~ "cmap" ]] ; then export BSSSI_CMAP=$FILE ; fi
                if [[ $FILE =~ "bnx" ]] ; then export BSSSI_BNX=$FILE ; fi
            fi
            if [[ $FILE =~ "DLE1" ]] ; then
                if [[ $FILE =~ "cmap" ]] ; then export DLE1_CMAP=$FILE ; fi
                if [[ $FILE =~ "bnx" ]] ; then export DLE1_BNX=$FILE ; fi
            fi
        done

        # report
        echo "BSPQI_CMAP: $BSPQI_CMAP"
        echo "BSPQI_BNX:  $BSPQI_BNX"
        echo "BSSSI_CMAP: $BSSSI_CMAP"
        echo "BSSSI_BNX:  $BSSSI_BNX"
        echo "DLE1_CMAP:  $DLE1_CMAP"
        echo "DLE1_BNX:   $DLE1I_BNX"

        # configuration for script identification
        export SUBMISSION_SCRIPT=""
        export CONFIGURATION_XML=""
        export ARGS=""

        # identify best data: DLE1+BssSi > DLE1 > BspQI+BssSI > BspQI  (currently DLE1+BssSI and BspQI are unsupported)
        if [[ ! -z "$DLE1_BNX" ]] && [[ ! -z "$DLE1_CMAP" ]] ; then
            if [[ ! -z "$BSSSI_BNX" ]] && [[ ! -z "$BSSSI_CMAP" ]] ; then
                echo "No support for DLE1 + BssSI, will just use DLE1."
            fi
            ln -s $DLE1_BNX DLE1
            ln -s $DLE1_CMAP DLE1.cmap
            export SUBMISSION_SCRIPT="hybrid_scaffold_3.3.sh"
            export CONFIGURATION_XML="hybridScaffold_DLE1_config_3.3.xml"
            export ARGS="DLE1"

        elif [[ ! -z "$BSPQI_BNX" ]] && [[ ! -z "$BSPQI_CMAP" ]] && [[ ! -z "$BSSSI_BNX" ]] && [[ ! -z "$BSSSI_CMAP" ]] ; then
            ln -s $BSPQI_BNX BSPQI
            ln -s $BSPQI_CMAP BSPQI.cmap
            ln -s $BSSSI_BNX BSSSI
            ln -s $BSSSI_CMAP BSSSI.cmap
            export SUBMISSION_SCRIPT="hybrid_scaffold_two.sh"
            export CONFIGURATION_XML="hybridScaffold_two_enzymes_largemem.xml"
            export ARGS="BSPQI BSSSI"

        else
            echo "Found no correct combination of input files!  Need both CMAP and BNX for either DLE1 or and BSPQI+BSSSI."
            exit 1
        fi


        # prepare to run script
        export OMP_NUM_THREADS=~{threadCount}
        export SLURM_CPUS_PER_TASK=~{threadCount}
        export ARGS="~{sampleName} $ARGS /root/config/$CONFIGURATION_XML"

        # run script
        bash /root/scripts/bionano/$SUBMISSION_SCRIPT $ARGS

	>>>
	output {
	    Array[File] outputFiles = glob(sampleName + "/*")
	}
    runtime {
        cpu: threadCount
        memory: "980 GB"
        docker: dockerRepository+"/vgp_bionano:"+dockerTag
    }
}