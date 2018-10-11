#! /bin/bash

if [[ -z $1 ]] || [[ ! -e input.fofn ]]; then
	echo "Usage: ./_submit_kat_comp.sh <asm.fasta>"
	echo "<asm.fasta>: Assembly fasta file"
	echo "<input.fofn>: Full path to the fastq files"
	exit -1
fi

asm=$1
input=input.fofn

for file in $(cat $input)
do
	# read-BC[1-2].fastq.gz
	file_name=`basename $file | sed 's/.*read/read/g'`
	if [ ! -e $file_name ]; then
		ln -s $file $file_name
	fi
done

asm_f=`basename $asm`

if [ ! -e $asm_f ]; then
	ln -s $asm 
fi
asm=$asm_f
genome=${asm/.fasta/}

cpus=32
mem=240g
name=$genome.kat
script=$VGP_PIPELINE/kat/kat_comp.sh
args="$asm"
partition=norm
walltime=2-0
path=`pwd`

mkdir -p logs
log=logs/$name.%A.log

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path --time=$walltime --error=$log --output=$log $script $args

