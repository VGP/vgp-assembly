#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: map.sh <asm.fasta> <input.fofn> [line_num]"
	echo "input.fofn: tab delimited file for R1 and R2 reads. ex. R1.fastq.gz\tR2.fastq.gz"
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

line=`sed -n ${i}p $input_fofn`

r1=`echo $line | awk '{print $1}'`
r2=`echo $line | awk '{print $2}'`

prefix="$i"

module load bwa/0.7.17
module load samtools

echo "
$tools/dfguan/utls/10x_trim -c -p $prefix $r1 $r2"
$tools/dfguan/utls/10x_trim -c -p $prefix $r1 $r2
echo

echo "align ${prefix}_1.fq.gz and ${prefix}_2.fq.gz to $prefix.bam"
echo "
bwa mem -t$SLURM_CPUS_PER_TASK $asm ${prefix}_1.fq.gz ${prefix}_2.fq.gz | samtools view -b - > $prefix.bam"
bwa mem -t$SLURM_CPUS_PER_TASK $asm ${prefix}_1.fq.gz ${prefix}_2.fq.gz | samtools view -b - > $prefix.bam
echo


if [ -s $prefix.bam ]; then
	echo "Found ${prefix}.bam."
	echo "
	rm ${prefix}_1.fq.gz ${prefix}_2.fq.gz"
	rm ${prefix}_1.fq.gz ${prefix}_2.fq.gz
fi
