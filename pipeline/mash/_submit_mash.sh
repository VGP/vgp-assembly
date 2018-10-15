#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_mash.sh <genome>"
	echo "Requires: input.fofn"
	exit -1
fi

LEN=`wc -l input.fofn | awk '{print $1}'`

genome=$1

mkdir -p logs

partition=quick
cpus=2
mem=2g
name=$genome.mash
script=$VGP_PIPELINE/mash/sketch.sh
args=""
walltime=240
path=`pwd`
log=logs/$name.%A_%a.log

echo "\
sbatch -J $name -a 1-$LEN --partition=$partition --cpus-per-task=$cpus -D $path --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name -a 1-$LEN --partition=$partition --cpus-per-task=$cpus -D $path --mem=$mem --time=$walltime --error=$log --output=$log $script $args > mash_jid

wait_for="--dependency=afterok:`cat mash_jid`"

cpus=4
mem=48g
name=$genome.plot
script=$VGP_PIPELINE/mash/plot.sh
walltime=30
args=$genome
log=logs/$name.%A_%a.log

echo "\
sbatch -J $name $wait_for --partition=$partition --cpus-per-task=$cpus -D $path --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name $wait_for --partition=$partition --cpus-per-task=$cpus -D $path --mem=$mem --time=$walltime --error=$log --output=$log $script $args

