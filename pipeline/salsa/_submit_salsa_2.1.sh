#!/bin/bash

if [[ -z $1 ]] || [[ ! -e re_bases.txt ]]  ; then
    echo "Usage: ./runArima.sh <fasta> [jobid]"
    echo -e "\tSymlink fastq.map : <R1.fastq.gz path> <R2.fastq.gz path>"
    echo -e "\tSymlink <fasta> : reference .fasta file to align"
    echo -e "\tSymlink re_bases.txt : Restriction enzyme site bases. ex. GATC,GANTC"
    exit -1
fi

ref=$1
ref_name=`basename $ref`
ref_name=`echo $ref_name | sed 's/.fasta$//g' | sed 's/.fa$//g'`

if ! [ -e $ref_name.fasta ]; then
	ln -s $ref $ref_name.fasta
fi
ref=$ref_name.fasta

path=`pwd`
cpus=4
mem=24g
name=$ref_name.index
log=logs/${name}.%A.log
jobid=$2
script=$VGP_PIPELINE/salsa/index.sh
partition=norm
extra=""
walltime=1-0
args=$ref

mkdir -p logs

echo "Run Salsa 2.1 against $ref"

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
fastq_map=fastq.map
log=logs/${name}.%A.log
script=$VGP_PIPELINE/salsa/arima_mapping_pipeline.sh
args="$fastq_map $ref_name $ref"
if ! [ -e $name.bed ]; then
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
mem=32g
name=$ref_name.salsa2.1
log=logs/${name}.%A.log
script=$VGP_PIPELINE/salsa/salsa2.1.sh
args="$ref_name ${ref_name}_salsa"
walltime=2-0

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args


