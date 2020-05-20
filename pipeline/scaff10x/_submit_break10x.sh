#!/bin/bash

if [[ -z $1 ]] || [[ ! -e asm.fasta ]] ; then
	echo "Usage: ./_submit_break10x.sh <genome_id> [job_id_to_wait_for]"
	echo "Requires asm.fasta, read-BC1.fastq and read-BC2.fastq"
	echo "[job_id_to_wait_for]: Wait after the job finishes and launch"
	exit -1
fi

mkdir -p logs

dependency=$2
if [ ! -z $dependency ]; then
	dependency="--dependency=afterok:$2"
fi

cpus=54
mem=64g
name=break10x_$1
script=$VGP_PIPELINE/scaff10x/break10x.sh
args=$1
walltime=3-0
dependency=$dependency
log=logs/$name.%A_%a.log
echo "sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
