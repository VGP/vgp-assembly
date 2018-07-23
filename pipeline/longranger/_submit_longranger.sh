#! /bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_longranger.sh <genome>"
	exit -1
fi

cpus=2
mem=12g
name=$1.longrgr
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/longranger/longranger.sh
args=$1
walltime=4-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

