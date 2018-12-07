#!/bin/bash

prefix=$1

if [ -z $SLURM_ARRAY_TASK_ID ]; then
	R=$2
else
	R=$SLURM_ARRAY_TASK_ID
fi

echo "\
cat $prefix.r$R.*.fastq | gzip > $prefix.read-BC$R.fastq.gz"
cat $prefix.r$R.*.fastq | gzip > $prefix.read-BC$R.fastq.gz
