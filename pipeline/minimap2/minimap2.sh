#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./minimap2 <ref>"
	echo "    Assumes we have <ref> and <input.fofn> in the same dir"
	exit -1
fi

module load minimap2/2.11
module load samtools

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

preset="$3"
if [ -z $3 ]; then
    preset="map-pb"
fi

samtools_filter="$4"
if [ -z $4 ]; then
    samtools_filter=""
fi

cpus=$SLURM_CPUS_PER_TASK

# Unless specified, use slurm array task id for input line num.
if [ -z $2 ]; then
	i=$SLURM_ARRAY_TASK_ID
else
	i=$2
fi

qry=`sed -n ${i}p input.fofn`

out=`basename $qry`
out=`echo $out | sed 's/.fasta.$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`
out=`echo $out | sed 's/.fastq.$//g' | sed 's/.fastq$//g' | sed 's/.fq$//g' | sed 's/.fastq.gz$//g' | sed 's/.fq.gz$//g'`


if [ -e $out.sort.bam ]; then
	echo "$out.sort.bam found. Exit."
	exit 0
fi

if [ -e $out.bam ]; then
	echo "$out.bam found. Skip alignment."
else
	echo "Start aligning $qry to $ref.idx"

	echo "\
	minimap2 -x $preset -a -t $cpus $ref.idx $qry | samtools view -hb $samtools_filter - > $out.bam"
	minimap2 -x $preset -a -t $cpus $ref.idx $qry | samtools view -hb $samtools_filter - > $out.bam
fi

echo "Sort $out.bam"

echo "\
samtools sort -T $out.tmp -O bam -o $out.sort.bam $out.bam"
samtools sort -T $out.tmp -O bam -o $out.sort.bam $out.bam
samtools index $out.sort.bam

rm $out.bam
