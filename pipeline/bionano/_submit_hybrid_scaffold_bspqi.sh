#!/bin/bash


if [ -z $1 ]; then
	echo "Usage: ./_submit_hybrid_scaffold_bspqi.sh <name>"
	echo -e "\t<name>: job name and output dir name"
	exit -1
fi

if [ -e $1 ]; then
	echo "$1 already exists. Remove."
	exit -1
fi

cpus=32
mem=180g
partition=largemem
name=$1
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/bionano/hybrid_scaffold.sh
args="$1 BSPQI /data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/bionano/hybridScaffold_config_BSPQI.xml"
walltime=3-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

