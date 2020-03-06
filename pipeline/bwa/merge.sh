#!/bin/bash

echo "Usage: ./merge.sh <out>"
echo "Merge bams in bam.list to <out>.bam"

cpus=$SLURM_CPUS_PER_TASK
out=$1

if [ -z $out ]; then
	exit -1
fi

module load samtools

bams=`cat bam.list | tr '\n' ' '`
bais=`echo $bams | sed 's/.bam/.bam.bai/g'`

num_bams=`wc -l bam.list | awk '{print $1}'`
if [[ "$num_bams" -eq 1 ]]; then
	echo "Only 1 bam provided. Skipping Merging."
	exit 0
fi

echo "Merge $bams"

echo "
samtools merge -@ $cpus -O bam -b bam.list $out.bam"
samtools merge -@ $cpus -O bam -b bam.list $out.bam

echo "
samtools index $out.bam"
samtools index $out.bam

echo "
rm $bams $bais"
