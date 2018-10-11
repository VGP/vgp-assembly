#!/bin/bash

if [[ -z $1 ]] || ! [[ -e input.fofn ]] ; then
	echo "Usage: ./_submit_juicer.sh <ref.fasta>"
	echo "Assumes input.fofn in the same path"
	echo "input.fofn should have reads with *_R1 and *_R2"
	exit -1
fi

ref=`realpath $1`

# Create reference dir
if ! [ -e reference/asm.fasta ]; then
	mkdir -p reference
	cd reference
	ln -s $ref asm.fasta
	cd ../
fi

# Create fastq dir
mkdir -p fastq
cd fastq
rm *
for fastq in $(cat ../input.fofn)
do
	ln -s $fastq
done
cd ../

ref=`basename $ref`
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g'`

cpus=4
mem=12g
name=juicer.$ref
script=$VGP_PIPELINE/juicer/juicer.sh
args=$ref
partition=norm
walltime=2-0
path=`pwd`
extra=$2

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args

