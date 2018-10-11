#!/bin/bash


if [ -z $1 ]; then
	echo "Usage: ./_submit_hybrid_scaffold_bspqi_bsssi.sh <name>"
	echo -e "\t<name>: job name and output dir name"
	exit -1
fi

cpus=64
mem=560g
partition=largemem
name=$1
script=$VGP_PIPELINE/bionano/hybrid_scaffold_two.sh
args="$1 BSPQI BSSSI $VGP_PIPELINE/bionano/hybridScaffold_two_enzymes_largemem.xml"
walltime=4-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

