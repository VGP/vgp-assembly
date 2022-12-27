#!/bin/bash

if [ -z $1 ]; then
  echo "Usage: ./_submit_telomere.sh <in.fasta>"
  echo
  echo "This script submits telomere_analysis.sh genome 0.4 50000 <in.fasta>"
  echo "    <threshold>: threshold to define telomere window. Define from 0.0 to 1.0. Recommended: 0.4"
  echo "    <ends>     : report ends.bed only if a window is found within <ends> bp of scaffold ends"
  exit -1
fi

fasta=$1
ln -s $fasta

fasta=`basename $fasta`
genome=`echo $fasta | sed 's/.fasta$//g' | sed 's/.fa$//g'`

cpus=12
mem=8g
name=$genome.telomere
#script=$VGP_PIPELINE/telomere/find_telomere.sh
script=$VGP_PIPELINE/telomere/telomere_analysis.sh
#args=$fasta
args="$genome 0.4 50000 $fasta"

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args

