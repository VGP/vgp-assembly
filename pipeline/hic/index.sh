#!/bin/sh

# Dependencies: bwa, samtools, picard

fa=$1

echo "bwa index $fa"
bwa index $fa -b 100000000

echo "samtools faidx $fa"
samtools faidx $fa

