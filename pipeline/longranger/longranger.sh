#!/bin/bash

if [ -z $1 ] ; then
	echo "Usage: ./longranger.sh <genome_id>"
	echo "Assumes we have asm.fasta"
	exit -1
fi

genome=$1
ref=asm.fasta
ref=${ref/.fasta/}      ## if contains .fasta, remove it

if ! [ -e "refdata-$ref/genome" ] ; then
	echo "No reference found"
	echo "=== start indexing reference ==="
	echo "\
	$tools/longranger/longranger-2.2.2/longranger mkref $ref.fasta"
	$tools/longranger/longranger-2.2.2/longranger mkref $ref.fasta
	echo ""
fi

echo "=== start running longranger wgs ==="
$tools/longranger/longranger-2.2.2/longranger align \
--id=$genome \
--fastq=/data/rhiea/genome10k/$genome/genomic_data/10x/ \
--sample=$genome \
--reference=refdata-$ref \
--jobmode=slurm \
--localcores=32 \
--localmem=60 \
--maxjobs=500 \
--jobinterval=5000 \
--disable-ui \
--nopreflight

