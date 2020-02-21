#!/bin/bash

ref=$1
fastq_map=$2
out=$3

if [[ -z $ref || -z $fastq_map || -z $out ]]; then
	echo "Usage: ./bwa.sh <ref.fasta> <fastq_map> <out>"
	echo "No <ref.fasta> found. Exit."
	exit -1
fi

module load bwa
module load samtools

cpus=$SLURM_CPUS_PER_TASK
idx=$SLURM_ARRAY_TASK_ID
out=$out.$idx

line=`sed -n ${idx}p $fastq_map`
r1=`echo $line | awk '{print $1}'`
r2=`echo $line | awk '{print $2}'`

echo "
bwa mem -t $cpus $ref $r1 $r2 > $out.sam"
bwa mem -t $cpus $ref $r1 $r2 > $out.sam

echo "
samtools sort -@$cpus -O bam -o $out.bam -T $out.tmp $out.sam"
samtools sort -@$cpus -O bam -o $out.bam -T $out.tmp $out.sam && rm $out.sam || exit -1

echo "
samtools index $out.bam"
samtools index $out.bam

