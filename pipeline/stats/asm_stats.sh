#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./asm_stats.sh <asm.fasta> <exp_genome_size (bp)>"
    exit -1
fi

fasta=$1
gsize=$2
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/stats/

asm=${fasta/.fasta/}

if ! [ -e $fasta.len ]; then
	java -jar -Xmx1g $script/fastaContigSize.jar $fasta
	java -jar -Xmx1g $script/lenCalcNGStats.jar $fasta.len $gsize > $asm.stats
fi
echo "Scaffolds"
cat $asm.stats

N_BASES=`awk '{sum+=$2; sumN+=$3} END {print (sum-sumN)}' $fasta.len`
echo "N bases: $N_BASES"
echo

if [ ! -e $asm.contigs.len ]; then
	java -jar -Xmx2g $script/fastaGetGaps.jar $fasta $asm.gaps
	awk -F "\t" '$4>3 {print $1"\t"$2"\t"$3}' $asm.gaps > $asm.gaps.bed
	awk '{print $1"\t0\t"$(NF-1)}' $fasta.len > $fasta.len.bed

	module load bedtools
	bedtools subtract -a $fasta.len.bed -b $asm.gaps.bed | awk '{print $NF-$(NF-1)}' > $asm.contigs.len
fi

echo "Contigs"
java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.contigs.len $gsize 1 > $asm.contigs.stats
cat $asm.contigs.stats

echo "Gaps"
java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.gaps $gsize 3 > $asm.gaps.stats
cat $asm.gaps.stats

rm $fasta.len.bed
rm $asm.gaps.bed
