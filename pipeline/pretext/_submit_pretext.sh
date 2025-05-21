#!/bin/bash

if [[ "$#" -lt 3 ]] ; then
    echo "Usage: ./_submit_pretext.sh <fasta> <fastq.map> <out_prefix> [jobid]"
	echo "  fasta     : reference .fasta file to align"
    echo "  fastq.map : /path/to/R1.fastq.gz <tab> /path/to/R2.fastq.gz"
	echo "  out_prefix: prefix for output files. Final output will be a merged <out_prefix>.bam and <out_prefix>.pretext"
    exit -1
fi

ref=$1

fastq_map=$2
if [ -z $fastq_map ]; then
	echo "No <fastq.map> provided. Exit."
	exit -1
fi

out_prefix=$3
if [ -z $out_prefix ]; then
	echo "No <out_prefix> provided. Exit."
	exit -1
fi

jobid=$4

path=`pwd`
cpus=4
mem=24g
name=$out_prefix.index
log=logs/${name}.%A.log
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
name=$out_prefix.map
walltime=1-0 # 3-0
gres="--gres=lscratch:100" # 600 for 600GB - for regular 50x HiC
log=logs/${name}.%A.log
script=$VGP_PIPELINE/salsa/arima_mapping_pipeline.sh
args="$fastq_map $out_prefix $ref"
if ! [ -e $out_prefix.bam ]; then
	echo "\
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra $gres --time=$walltime --error=$log --output=$log $script $args"
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra $gres --time=$walltime --error=$log --output=$log $script $args >> mapping_jid
	jobid=`tail -n1 mapping_jid`

fi

if ! [ -z $jobid ] ; then
	extra="--dependency=afterok:$jobid"
else
	extra=""
fi

cpus=4
mem=4g
name=$out_prefix.pretext
log=logs/${name}.%A.log
script=$VGP_PIPELINE/pretext/pretext.sh
args="$out_prefix.bam"
walltime=1-0

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args


