#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./_submit_freebayes.sh <sample_id> [jobid_to_set_dependency]"
    exit -1
fi


sample=$1

if ! [ -e aligned.bam ]; then
	ln -s $sample/outs/possorted_bam.bam aligned.bam
	ln -s $sample/outs/possorted_bam.bam.bai aligned.bam.bai
fi

mkdir -p logs
cpus=2
mem=8g
name=$1.freebayes
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/freebayes-polish/freebayes.sh
args=$sample
walltime=8:00:00
log=logs/$name.%A_%a.log

mkdir -p bcf

if ! [ -z $2 ]; then
	wait_for="--dependency=afterok:$2"
else
	echo "\
	sbatch --partition=norm --array=1-100 -D $PWD --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=norm --array=1-100 -D $PWD --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > freebayes_jid
	wait_for="--dependency=afterok:`cat freebayes_jid`"
fi

cpus=2
mem=4g
name=$1.consensus
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/freebayes-polish/consensus.sh
args=$sample
walltime=2-0
log=logs/$name.%A_%a.log

echo "\
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

