#! /bin/bash

echo "Usage: ./_submit_bwa.sh <ref.fasta> <fastq.map> <out>"

ref=$1
fastq_map=$2
out=$3

if [ -z $ref ]; then
	echo "No <ref.fasta> given. Exit."
	exit -1
fi

mkdir -p logs

if [ ! -s $ref.fai ]; then
	cpus=4
	mem=10g
	name=$out.index
	script=$VGP_PIPELINE/bwa/bwa_index.sh
	args="$ref"
	partition=quick
	walltime=4:00:00
	path=`pwd`
	log=logs/$name.%A.log

	echo "\
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > index.jid
	index=`cat index.jid`
	echo $index
	extra="--dependency=$index"
fi

cpus=24
mem=120g
name=$out.bwa
script=$VGP_PIPELINE/bwa/bwa.sh
args="$ref $fastq_map $out"
partition=norm
walltime=1-0
path=`pwd`
log=logs/$name.%A_%a.log

LEN=`wc -l $fastq_map | awk '{print $1}'`
extra="--array=1-$LEN $extra"

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > bwa.jid
wait_for="--dependency=afterok:"`cat bwa.jid`

for i in $(seq 1 $LEN)
do
	echo "$out.$i.bam" >> bam.list
done

cpus=24
mem=48g
name=$out.merge
script=$VGP_PIPELINE/bwa/merge.sh
args="$out"
partition=norm
walltime=1-0
path=`pwd`
log=logs/$name.%A.log
extra=$wait_for

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args

