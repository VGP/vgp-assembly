#!/bin/bash

if [[ -z $1 ]] ; then
	echo "Usage: ./_submit_hic.sh <name>"
	echo "	Run with an existing .bam file"
	echo "  Requires gaps.bed"
	exit -1
fi

#export asset=$tools/asset
name=$1

mkdir -p logs

if [[ ! -e gaps.bed ]]; then
	echo "No gaps.bed found. Exit."
	exit -1
fi

if [[ ! -e hic.bed ]]; then
	cpus=4
	mem=10g
	name=hic.ast.$name
	script=$asset/scripts/slurm/hic/ast.sh
	args=gaps.bed
	walltime=4:00:00
	dependency=""
	log=$PWD/logs/$name.%A.log
	echo "\
	sbatch --partition=quick -D `pwd` --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=quick -D `pwd` --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
else
	echo "*** hic.bed found. nothing will be run. ***"
fi

