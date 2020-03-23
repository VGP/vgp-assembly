#!/bin/bash

echo "Usage: ./pretext.sh <in.bam>"
echo "  This will generate <in>.pretext"

module load samtools

bam=$1

if [ -z $bam ]; then
	echo "No bam provided. Exit."
	exit -1
fi

out=${bam/.bam/.pretext}

echo "\
samtools view -h $bam | $tools/Pretext/PretextMap -o $out"
samtools view -h $bam | $tools/Pretext/PretextMap -o $out
