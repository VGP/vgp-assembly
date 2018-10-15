#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_busco.sh <assembly.fasta>"
	exit -1
fi

assembly=$1
if ! [ -e `basename $assembly` ]; then
	ln -s $assembly
fi

assembly=`basename $assembly`
name=${assembly/.fasta/}

cpus=24
mem=42g
name=$name.busco
script=$VGP_PIPELINE/busco/busco.sh
args=$assembly
partition=norm
walltime=4-0
path=`pwd`

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --time=$walltime --error=$log --output=$log $script $args

