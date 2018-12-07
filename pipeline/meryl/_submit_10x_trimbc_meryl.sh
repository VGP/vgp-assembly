#!/bin/bash

if [ -z $1 ]; then
        echo "Usage: sh _submit_10x_trimbc_meryl.sh <k-size> <genome_id>"
	echo "requires input.fofn: ls *R1*.fastq.gz > input.fofn"
        exit -1
fi


k=$1	# 31
genome_id=$2	# genome_id

if ! [ -e input.fofn ]; then
	echo "requires input.fofn: ls --color=never *R1_001.fastq* > input.fofn"
	exit -1
fi

mkdir -p logs

## Re-header 10X fastq files with scaff10x

LEN=`wc -l input.fofn | awk '{print $1}'`
path=$PWD

partition=norm
cpus=4
mem=40g
name=trim_bc.$genome_id
log=logs/$name.%A_%a.log
walltime=2-0
script=$VGP_PIPELINE/meryl/trim_bc.sh
args="$genome_id"

echo "\
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime --error=$log --output=$log -D $path --array=1-$LEN $script $args"
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime --error=$log --output=$log -D $path --array=1-$LEN $script $args > trim_bc.jid

jid=`cat trim_bc.jid`
if [ -z $jid ]; then
	dependency=""
else
	dependency="--dependency=afterok:$jid"
fi

if ! [ -e fastq.list ]; then
	for i in $(seq 1 $LEN); do
		echo "$genome_id.r1.$i.fastq" >> fastq.list
		echo "$genome_id.r2.$i.fastq" >> fastq.list
	done
fi

## Concatenate re-headered fastq files

partition=norm
cpus=4
mem=8g
walltime="2-0"
array=1-2
name=concat.$genome_id
script=$VGP_PIPELINE/meryl/concat.sh
log=logs/$name.%A_%a.log
args="$genome_id"

echo "\
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime $dependency --error=$log --output=$log -D $path --array=$array $script $args"
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime $dependency --error=$log --output=$log -D $path --array=$array $script $args

LEN=`wc -l fastq.list | awk '{print $1}'`
partition=largemem
cpus=24
mem=120g
walltime="2-0"
name=meryl_build.$genome_id
log=logs/$name.%A_%a.log
script=$VGP_PIPELINE/meryl/meryl_build.sh
args="$k fastq.list"

echo "\
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime $dependency --error=$log --output=$log -D $path --array=1-$LEN $script $args"
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime $dependency --error=$log --output=$log -D $path --array=1-$LEN $script $args > meryl_build.jid
jid=`cat meryl_build.jid`
if [ -z $jid ]; then
        dependency=""
else
        dependency="--dependency=afterok:$jid"
fi

cpus=4
mem=60g
partition=norm
walltime="1-0"
name=meryl_mrg.$genome_id
log=logs/${name}_%A.log
script=$VGP_PIPELINE/meryl/meryl_mrg.sh
args="$k $genome_id fastq.list"

echo "\
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime $dependency --error=$log --output=$log -D $path $script $args"
sbatch -J $name -c $cpus --mem=$mem --time=$walltime -P $partition -t $walltime $dependency --error=$log --output=$log -D $path $script $args > meryl_mrg.jid
