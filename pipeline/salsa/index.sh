#!/bin/sh

module load bwa/0.7.17
module load samtools/1.9

fa=$1

echo "bwa index $fa"
bwa index $fa
echo

echo "samtools faidx $fa"
samtools faidx $fa
echo

