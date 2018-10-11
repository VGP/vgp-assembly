#! /bin/bash

if [[ -z $1 ]] ; then
	echo "Usage: ./_submit_indexing.sh <genome> <ref.fasta>"
	echo "<ref.fasta> will be linked to asm.fasta."
	echo "This script is doing indexing only, in case it failed for large genomes."
	exit -1
fi

if ! [ -e asm.fasta ]; then
	ln -s $2 asm.fasta
fi

cpus=8
mem=32g
name=$1.longrgr
script=$VGP_PIPELINE/longranger/indexing.sh
args=$1
walltime=2-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

