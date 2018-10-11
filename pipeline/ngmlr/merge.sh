#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: _mrg.sh <out_prefix>"
	exit -1
fi

module load samtools
prefix=$1
echo "\
samtools merge -O bam -@$SLURM_CPUS_PER_TASK $prefix.bam *.sorted.bam"
samtools merge -O bam -@$SLURM_CPUS_PER_TASK $prefix.bam *.sorted.bam
echo "\
samtools index $prefix.bam"
samtools index $prefix.bam
