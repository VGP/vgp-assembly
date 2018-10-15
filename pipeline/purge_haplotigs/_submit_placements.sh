#! /bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_placements.sh <genome>"
	exit -1
fi

genome=$1

ln -s ../../${genome}_v1.p.fasta
ln -s ../../${genome}_v1.h.fasta

cpus=24
mem=48g
name=$genome.placements
script=$VGP_PIPELINE/purge_haplotigs/placements.sh
args=$genome
partition=norm
walltime=2-0
path=`pwd`

mkdir -p logs
log=logs/$name.%A_%a.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --time=$walltime --error=$log --output=$log $script $args

