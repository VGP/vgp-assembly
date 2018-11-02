#!/bin/bash

if [[ -z $1 ]] || [[ ! -e re_bases.txt ]]  ; then
    echo "Usage: ./runArima.sh <name> [jobid]"
    echo -e "\tSymlink fastq.map : <R1.fastq.gz path> <R2.fastq.gz path>"
    echo -e "\tSymlink <name>.fasta : reference fasta file to align"
    echo -e "\tSymlink re_bases.txt : Restriction enzyme site bases. ex. GATC,GANTC"
    exit -1
fi

threads=48
mem=40g
name=$1
fastq_map=fastq.map
ref=$name.fasta
log=logs/map_${name}.%A.%a.log
jobid=$2
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/salsa/arima_mapping_pipeline.sh

bwa_indexing=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/salsa/index.sh

wait_until=""

mkdir -p logs

if ! [ -z $jobid ] ; then
	    wait_until="--dependency=afterok:$jobid"
elif [ ! -e $ref.fai ]; then
	echo "\
	sbatch --partition=quick --time=4:00:00 --cpus-per-task=2 --mem=4g --error=logs/index_${name}.%A.log --output=logs/_${name}.%A.log --job-name=$name.index $bwa_indexing $ref > logs/index_$name.log"
	sbatch --partition=quick --time=4:00:00 --cpus-per-task=2 --mem=4g --error=logs/index_${name}.%A.log --output=logs/_${name}.%A.log --job-name=$name.index $bwa_indexing $ref > logs/index_$name.log
	jobid=`cat logs/index_$name.log | tail -n 1`
	wait_until="--dependency=afterok:$jobid"
fi

if ! [ -e $name.bed ]; then
	echo "\
	sbatch --partition=norm --time=2-0 --cpus-per-task=$threads --mem=$mem --gres=lscratch:500 --error=$log --output=$log $wait_until --job-name=$name $script $fastq_map $name $ref $threads"
	sbatch --partition=norm --time=2-0 --cpus-per-task=$threads --mem=$mem --gres=lscratch:500 --error=$log --output=$log $wait_until --job-name=$name $script $fastq_map $name $ref $threads > mapping_jid
	jobid=`cat mapping_jid`

fi

if ! [ -z $jobid ] ; then
	wait_until="--dependency=afterok:$jobid"
else
	wait_until=""
fi

threads=2
mem=32g
name=${1}
log=logs/salsa_${name}.%A.%a.log
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/salsa/salsa2.sh

echo "\
sbatch --partition=norm --time=2-0 --cpus-per-task=$threads --mem=$mem --error=$log --output=$log $wait_until --job-name=$name $script $1 ${name}_salsa"
sbatch --partition=norm --time=2-0 --cpus-per-task=$threads --mem=$mem --error=$log --output=$log $wait_until --job-name=$name $script $1 ${name}_salsa


