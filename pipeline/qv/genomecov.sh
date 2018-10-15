#!/bin/bash

module load bedtools
module load samtools

if [ -z $1 ]; then
	echo "Usage: ./genomecov.sh <genome_id> [bam]"
	exit -1
fi

bam=aligned.bam		#aligned, pos-sorted
genome=$1

if [ ! -z $2 ]; then
	bam=$2
fi

if [ -e aligned.genomecov ]; then
	echo "aligned.genomecov already exists. skip..."
else
	echo "Collect coverage"
	echo "\
	samtools view -F 0x100 -hb $bam | bedtools genomecov -ibam - -split > aligned.genomecov"
	samtools view -F 0x100 -hb $bam | bedtools genomecov -ibam - -split > aligned.genomecov
fi
echo

echo "\
awk '{if ($1=="genome" && $2>3) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp"
awk '{if ($1=="genome" && $2>3) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp
NUM_BP=`cat $genome.numbp`
echo "Total bases > 3x: echo $NUM_BP"

if [ -e $genome.numvar ]; then
	NUM_VAR=`cat $genome.numvar`
	echo "Total num. bases subject to change: $NUM_VAR"
	QV=`echo "$NUM_VAR $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}'`
	echo $QV > $genome.qv
	echo "QV of this genome $genome: $QV"
fi
