#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./_submit_nucmer.sh <ref.fasta> <qry.fasta>"
    exit 0
fi

ref=$1
qry=$2

ref_prefix=`basename $ref`
ref_prefix=`echo $ref_prefix | sed 's/.fasta$//g' | sed 's/.fa$//g'`
qry_prefix=`basename $qry`
qry_prefix=`echo $qry_prefix | sed 's/.fasta$//g' | sed 's/.fa$//g'`

cpus=12
mem=100g
name=nucmer.$qry_prefix.to.$ref_prefix
script=$VGP_PIPELINE/nucmer/nucmer.sh
args="$ref $qry"
partition=norm
walltime=1-0
path=`pwd`
log=logs/$name.%A.log

mkdir -p logs

echo "\
sbatch -J $name --mem=$mem --cpus-per-task=$cpus -D $path --partition=$partition --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --cpus-per-task=$cpus -D $path --partition=$partition --time=$walltime --error=$log --output=$log $script $args

