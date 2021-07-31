#!/bin/bash

if [[ -z $1 ]] ; then
    echo "Usage: ./_submit_pretext.sh <fasta> <fastq.map> [jobid]"
    echo -e "\t<fastq.map> : <R1.fastq.gz path> <R2.fastq.gz path>"
    echo -e "\tSymlink <fasta> : reference .fasta file to align"
    exit -1
fi

ref=$1
ref_name=`basename $ref`
ref_name=`echo $ref_name | sed 's/.fasta$//g' | sed 's/.fa$//g'`

fastq_map=$2
if [ -z $fastq_map ]; then
	echo "No <fastq.map> provided. Exit."
	exit -1
fi

if ! [ -e $ref_name.fasta ]; then
	ln -s $ref $ref_name.fasta
fi
ref=$ref_name.fasta

path=`pwd`
cpus=4
mem=24g
name=$ref_name.index
log=logs/${name}.%A.log
jobid=$3
script=$VGP_PIPELINE/salsa/index.sh
partition=norm
extra=""
walltime=1-0
args=$ref

mkdir -p logs

echo "Run Arima mapping pipeline against $ref"

if ! [ -z $jobid ] ; then
	    extra="--dependency=afterok:$jobid"
elif [ ! -e $ref.fai ]; then
	echo "\
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > index_jid
	jobid=`cat index_jid | tail -n 1`
	extra="--dependency=afterok:$jobid"
fi

cpus=48
mem=40g
name=$ref_name.map
walltime=3-0
log=logs/${name}.%A.log
script=$VGP_PIPELINE/salsa/arima_mapping_pipeline.sh
args="$fastq_map $ref_name $ref $cpus"
if ! [ -e $ref_name.bam ]; then
	echo "\
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --gres=lscratch:600 --time=$walltime --error=$log --output=$log $script $args"
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --gres=lscratch:600 --time=$walltime --error=$log --output=$log $script $args > mapping_jid
	jobid=`cat mapping_jid`

fi

if ! [ -z $jobid ] ; then
	extra="--dependency=afterok:$jobid"
else
	extra=""
fi

cpus=4
mem=4g
name=$ref_name.pretext
log=logs/${name}.%A.log
script=$VGP_PIPELINE/pretext/pretext.sh
args="$ref_name.bam"
walltime=1-0

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args


