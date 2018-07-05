#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./run.sh <genome_id>"
	exit -1
fi

mkdir -p logs

cpus=48
mem=48g
name=scaff10x_$1
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/scaff10x/scaff10x.sh
args=$1
walltime=3-0
dependency=""
log=logs/$name.%A_%a.log
echo "sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
