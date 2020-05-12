version 1.0

workflow runBionanoSolve {
	call bionano_solve
}

task bionano_solve {
    input {
        Array[File] bionanoFiles
        File assemblyFasta
        String sampleName
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

        # potential files and cmaps
        export BSPQI_CMAP=""
        export BSSSI_CMAP=""
        export DLE1_CMAP=""

        # look through all files
        shopt -s nocasematch
        for FILE in ~{sep=" " bionanoFiles} ; do
            if [[ $FILE =~ "BspQI" ]] && [[ $FILE =~ "cmap" ]] ; then
                export BSPQI_CMAP=$FILE
            fi
            if [[ $FILE =~ "BssSI" ]] && [[ $FILE =~ "cmap" ]] ; then
                export BSSSI_CMAP=$FILE
            fi
            if [[ $FILE =~ "DLE1" ]] && [[ $FILE =~ "cmap" ]] ; then
                export DLE1_CMAP=$FILE
            fi
        done

        # report
        echo "BSPQI_CMAP: $BSPQI_CMAP"
        echo "BSSSI_CMAP: $BSSSI_CMAP"
        echo "DLE1_CMAP:  $DLE1_CMAP"

        # configuration for script identification
        export SUBMISSION_SCRIPT=""
        export CONFIGURATION_XML=""
        export ARGS=""

        # identify best data: DLE1+BssSi > DLE1 > BspQI+BssSI > BspQI  (currently DLE1+BssSI and BspQI are unsupported)
        if [[ ! -z "$DLE1_CMAP" ]] ; then
            if [[ ! -z "$BSSSI_CMAP" ]] ; then
                echo "No support for DLE1 + BssSI, will just use DLE1."
            fi
            ln -s $DLE1_CMAP DLE1.cmap
            export SUBMISSION_SCRIPT="hybrid_scaffold_3.3.sh"
            export CONFIGURATION_XML="hybridScaffold_DLE1_config_3.3.xml"
            export ARGS="DLE1"

        elif [[ ! -z "$BSPQI_CMAP" ]] && [[ ! -z "$BSSSI_CMAP" ]] ; then
            ln -s $BSPQI_CMAP BSPQI.cmap
            ln -s $BSSSI_CMAP BSSSI.cmap
            export SUBMISSION_SCRIPT="hybrid_scaffold_two.sh"
            export CONFIGURATION_XML="hybridScaffold_two_enzymes_largemem.xml"
            export ARGS="BSPQI BSSSI"

        else
            echo "Found no correct combination of input files!  Need CMAP for DLE1 or BSPQI+BSSSI."
            exit 1
        fi


        # prepare to run script
        export OMP_NUM_THREADS=~{threadCount}
        export SLURM_CPUS_PER_TASK=~{threadCount}
        export ARGS="~{sampleName} $ARGS /root/config/bionano/ $CONFIGURATION_XML"
        export tools="/opt"
        ln -s ~{assemblyFasta} asm.fasta

        # run script bionano script
        bash /root/scripts/bionano/$SUBMISSION_SCRIPT $ARGS

        # get assembly data
        cd ~{sampleName}
        tar xvf hybridScaffold_archive.tar.gz
        cat hybridScaffold_archive/hybrid_scaffolds/*_HYBRID_SCAFFOLD_NCBI.fasta hybridScaffold_archive/hybrid_scaffolds/*_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta > ../scaffolded_asm.fasta
        cd ..

        # finalize and trim N's
        bash /root/scripts/bionano/trimNs/trimNs.sh

	>>>
	output {
	    File bionanoTarball = sampleName + "/hybridScaffold_archive.tar.gz"
	}
    runtime {
        cpu: threadCount
        memory: "980 GB"
        docker: dockerImage
    }
}