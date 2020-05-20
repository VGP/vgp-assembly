#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: map_hic.sh <asm.fasta> <input.fofn> [line_num num_cpus]"
	echo -e "\tinput.fofn: tab delimited file for R1 and R2 reads. ex. R1.fastq.gz\tR2.fastq.gz"
	echo -e "\t[line_num]: line num. to proceed in <input.fofn>"
	echo -e "\t[num_cpus]: num. cpus used (if not running on SLURM)"
	exit -1
fi

asm=$1
input_fofn=$2

i=$SLURM_ARRAY_TASK_ID
if [ -z $i ]; then
	i=$3
fi

if [ -z $i ]; then
	echo "No line num provided. Exit."
	exit -1
fi

cpus=$SLURM_CPUS_PER_TASK
if [ -z $cpus ]; then
	cpus=$4
fi

line=`sed -n ${i}p $input_fofn`

r1=`echo $line | awk '{print $1}'`
r2=`echo $line | awk '{print $2}'`

prefix="$i"

module load bwa/0.7.17
module load samtools

echo "align $r1.fq.gz and $r2.fq.gz to $prefix.bam"
echo "
bwa mem -SP -B10 -t$cpus $asm $r1 $r2 | samtools view -b - > $prefix.bam"
bwa mem -SP -B10 -t$cpus $asm $r1 $r2 | samtools view -b - > $prefix.bam
echo


