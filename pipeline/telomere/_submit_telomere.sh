#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_telomere.sh <in.fasta>"
	exit -1
fi

fasta=$1
ln -s $fasta

fasta=`basename $fasta`
fasta_prefix=`echo $fasta | sed 's/.fasta$//g' | sed 's/.fa$//g'`

cpus=12
mem=8g
name=$fasta_prefix.telomere
script=$VGP_PIPELINE/telomere/find_telomere.sh
args=$fasta

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args

