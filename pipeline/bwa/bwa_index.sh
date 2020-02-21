#!/bin/bash

ref=$1

if [ -z $ref ]; then
	echo "Usage: ./bwa_index.sh <ref.fasta>"
	echo "No <ref.fasta> given. Exit."
	exit -1
fi

module load bwa
module load samtools

echo "
bwa index $ref"
bwa index $ref

echo "
samtools faidx $ref"
samtools faidx $ref

