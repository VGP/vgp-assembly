#!/bin/bash


if [ -z $1 ]; then
	echo "Usage: ./_submit_cmap_assembly.sh <name> <bnx> <out> <xml>"
	echo -e "\t<name>: job name"
	echo -e "\t<bnx>: path to the bnx file"
	echo -e "\t<out>: output folder for this bnx cmap assembly"
	echo -e "\t<xml>: /bionano_path/RefAligner/7437.7523rel/optArguments_<haplotype>_<instrument>....xml"
	exit -1
fi

cpus=52
mem=120g
partition=norm
name=$1
script=$VGP_PIPELINE/bionano/cmap_assembly_ref.sh
args="$2 $3 $4"
walltime=10-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=$partition --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

