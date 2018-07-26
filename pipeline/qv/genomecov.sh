#!/bin/bash

module load bedtools
module load samtools


bam=aligned.bam		#aligned, pos-sorted
genome=$1

echo "Collect coverage"
echo "\
samtools view -F 0x100 -hb $bam | bedtools genomecov -ibam - -split > aligned.genomecov"
samtools view -F 0x100 -hb $bam | bedtools genomecov -ibam - -split > aligned.genomecov
echo

echo "\
awk '{if ($1=="genome" && $2>3) {numbp += $3}} END {print numbp}' > $genome.numbp"
awk '{if ($1=="genome" && $2>3) {numbp += $3}} END {print numbp}' > $genome.numbp

NUM_BP=`cat $genome.numbp`
NUM_VAR=`cat $genome.numvar`
echo "Total bases > 3x: echo $NUM_BP"
echo "Total num. bases subject to change: $NUM_VAR"
QV=`echo "$NUM_VAR $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}'`
echo $QV > $genome.qv
echo "QV of this genome $genome: $QV"

