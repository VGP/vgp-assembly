#! /bin/bash

if [[ -z $1 ]] ; then
	echo "Usage: ./_submit_repeatmasker.sh <asm.fasta>"
	exit -1
fi

if ! [ -e asm.fasta ]; then
	ln -s $1 asm.fasta
fi

asm=`basename $1`
asm=${asm/.fasta/}

echo "Run repeatmasker on $asm using species=vertebrates"

cpus=54
mem=140g
partition=unlimited
name=$asm.rm
script=$VGP_PIPELINE/repeatmasker/repeatmasker.sh
args=""
walltime=10-0
path=`pwd`

mkdir -p logs
log=logs/$name.%A_%a.log

echo "\
sbatch --partition=$partition --cpus-per-task=$cpus -D $path -J $name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=$partition --cpus-per-task=$cpus -D $path -J $name --mem=$mem --time=$walltime --error=$log --output=$log $script $args

