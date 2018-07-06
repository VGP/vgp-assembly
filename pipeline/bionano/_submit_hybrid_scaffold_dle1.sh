#!/bin/bash


if [ -z $1 ]; then
	echo "Usage: ./_submit_hybrid_scaffold_dle1.sh <name>"
	echo -e "\t<name>: job name and output dir name"
	exit -1
fi

if [ -e $1 ]; then
	echo "$1 already exists. Remove."
	exit -1
fi

cpus=8
mem=80g
partition=norm
name=$1
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/bionano/hybrid_scaffold.sh
args="$1 DLE1 /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/hybridScaffold_DLE1_config.xml"
walltime=4-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

