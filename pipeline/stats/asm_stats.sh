#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./asm_stats.sh <asm.fasta> <exp_genome_size (bp)> [p/c]"
    echo "[p]: set for getting primary / alt haplotig stats, assuming the primary begins with scaffold_ in its name"
    echo "[c]: set for getting scaffolds (direct sttas) only."
    exit -1
fi

fasta=$1
gsize=$2
script=$VGP_PIPELINE/stats/

asm=${fasta/.fasta/}

if ! [ -e $fasta.len ]; then
	java -jar -Xmx1g $script/fastaContigSize.jar $fasta
fi
echo "Scaffolds"
java -jar -Xmx1g $script/lenCalcNGStats.jar $fasta.len $gsize > $asm.stats
cat $asm.stats

c=$3
if [[ "$c" == "c" ]]; then
        exit 0
fi

N_BASES=`awk '{sum+=$2; sumN+=$3} END {print (sum-sumN)}' $fasta.len`
echo "N bases: $N_BASES"
echo

if [ ! -e $asm.contigs.len ]; then
	java -jar -Xmx2g $script/fastaGetGaps.jar $fasta $asm.gaps
	awk -F "\t" '$4>3 {print $1"\t"$2"\t"$3}' $asm.gaps > $asm.gaps.bed
	awk '{print $1"\t0\t"$(NF-1)}' $fasta.len > $fasta.len.bed

	module load bedtools
	bedtools subtract -a $fasta.len.bed -b $asm.gaps.bed | awk '{print $1"\t"$NF-$(NF-1)}' > $asm.contigs.len
fi

echo "Contigs"
#if [ ! -e $asm.contigs.stats ]; then
	java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.contigs.len $gsize 1 > $asm.contigs.stats
#fi
cat $asm.contigs.stats

echo "Gaps"
#if [ ! -e $asm.gaps.stats ]; then
	java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.gaps $gsize 3 > $asm.gaps.stats
#fi
cat $asm.gaps.stats

#rm $fasta.len.bed
#rm $asm.gaps.bed

p=$3
if [[ "$p" != "p" ]]; then
	exit 0
fi
echo

echo "=== Primary Stats ==="
grep "scaffold_" $fasta.len > $asm.p.len
java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.p.len $gsize > $asm.p.stats

echo "Scaffolds"
cat $asm.p.stats
echo

grep "scaffold_" $asm.contigs.len > $asm.contigs.p.len
java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.contigs.p.len $gsize 1 > $asm.contigs.p.stats
echo "Contigs"
cat $asm.contigs.p.stats
echo

grep "scaffold_" $asm.gaps > $asm.gaps.p
java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.gaps.p $gsize 3 > $asm.gaps.p.stats
echo "Gaps"
cat $asm.gaps.p.stats
echo

if [ ! -e $asm.p.fasta ]; then
    echo "Extract primary set"
    cut -f1 $asm.p.len > $asm.p.list
    java -jar -Xmx1g $script/fastaExtractFromList.jar $asm.fasta $asm.p.list $asm.p.fasta
fi
echo

echo "=== Alt Stats ==="
grep -v "scaffold_" $fasta.len > $asm.h.len
java -jar -Xmx1g $script/lenCalcNGStats.jar $asm.h.len $gsize > $asm.h.stats
cat $asm.h.stats
echo

#:<<'END'
if [ ! -e $asm.h.fasta ]; then
    echo "Extract alt set"
    cut -f1 $asm.h.len > $asm.h.list
    java -jar -Xmx1g $script/fastaExtractFromList.jar $asm.fasta $asm.h.list $asm.h.fasta
fi
#END
