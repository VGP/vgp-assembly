version 1.0

workflow runSalsa {
	call salsa
}

task salsa {
    input {
        File refFasta
        Array[File] readFilesHiC
        String sampleName
        Int threadCount
        String? enzymeBases="GATC,GANTC"
        String? bwaMemOpts="-B8"
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

        # input reference
        REF=`basename ~{refFasta}`
        ln -s ~{refFasta}

        # prepare input read files
        touch R1.fastq.gz
        touch R2.fastq.gz
        for RF in ~{sep=" " readFilesHiC} ; do
            # should not happen, but just in case
            if [[ ! $RF == *gz ]] ; then
                echo "HiC files expected to be in .gz format, zipping: $RF"
                gzip -k $RF
                RF="$RF.gz"
            fi

            # sort into R1 and R2
            if [[ $RF == *R1* ]] && [[ $RF == *R2* ]] ; then
                echo "Cannot determine R1 or R2 from filename (found both): $RF"
                exit 1
            elif [[ $RF == *R1* ]] ; then
                cat $RF >> R1.fastq.gz
            elif [[ $RF == *R2* ]] ; then
                cat $RF >> R2.fastq.gz
            else
                echo "Cannot determine R1 or R2 from filename (found neither): $RF"
                exit 1
            fi
        done
        if [[ ! -s R1.fastq.gz ]] || [[ ! -s R2.fastq.gz ]] ; then
            echo "Empty R1 and R2 files"
            exit -1
        fi
        echo "R1.fastq.gz R2.fastq.gz" >>fastq.map

        # run preparation scripts
        export SLURM_CPUS_PER_TASK=~{threadCount}
        export VGP_PIPELINE="/root/scripts"
        export SLURM_JOBID="tmp"
        bash /root/scripts/salsa/index.sh $REF
        bash /root/scripts/salsa/arima_mapping_pipeline.sh fastq.map ~{sampleName} $REF `pwd` "~{bwaMemOpts}"
        module load bedtools
        bedtools bamtobed -i ~{sampleName}.bam > ~{sampleName}.bed

        # run salsa pipeline
        module load python/2.7
        mkdir out
        python /root/tools/salsa/SALSA-2.2/run_pipeline.py -a $REF -l $REF.fai -e ~{enzymeBases} -b ~{sampleName}.bed -o out -m yes -p yes
        mv out/scaffolds_FINAL.fasta ~{sampleName}.salsa.fasta

        # cleanup
        rm -r tmp R1.fastq.gz R2.fastq.gz

	>>>
	output {
		File scaffoldedFasta = sampleName + ".salsa.fasta"
	}
    runtime {
        cpu: threadCount
        memory: "40 GB"
        docker: dockerImage
#        docker: dockerRepository+"/vgp_salsa:"+dockerTag
    }
}
