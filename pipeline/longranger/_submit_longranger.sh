#! /bin/bash

if [[ -z $1 ]] ; then
	echo "Usage: ./_submit_longranger.sh <genome> <ref.fasta>"
	echo "<ref.fasta> will be linked to asm.fasta."
	echo "Assumes we have 10x reads located under /data/rhiea/<genome>/genomic_data/10x/."
	exit -1
fi

if ! [ -e asm.fasta ]; then
	ln -s $2 asm.fasta
fi

cpus=2
mem=12g
name=$1.longrgr
script=$VGP_PIPELINE/longranger/longranger.sh
args=$1
walltime=4-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

