#!/bin/bash

if [[ -z $1 ]] ; then
	echo "Usage: ./_submit_10x.sh <genome> [mean_cov]"
	echo "Asset will run on aligned.bam from longranger and gaps.bed in this dir."
	exit 0;
fi

if [ ! -e aligned.bam ]; then
	echo "No \"aligned.bam\" found."
	exit -1;
fi

genome=$1
mean_cov=$2

mkdir -p logs

#if [ ! -e 10x.bed ]; then
	cpus=4
	mem=48g
	name=10x.ast.$genome
	script=$VGP_PIPELINE/asset/longranger/ast.sh
	args="gaps.bed $mean_cov"
	walltime=4:00:00
	log=$PWD/logs/$name.%A.log
	echo "\
	sbatch --partition=quick --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=quick --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
#fi

