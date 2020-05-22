#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./asm_stats.sh <asm.fasta> <exp_genome_size (bp)> [p/s/c]"
    echo "[ ]: By default, stats for scaffold / contig / gaps will be generated."
    echo "[p]: set for getting primary / alt haplotig stats, in addition to default, assuming the primary begins with scaffold_ in its name"
    echo "[s]: set for getting scaffold stats (direct stats) only."
    echo "[c]: set for getting contig (direct stats) as well."
    exit -1
fi

#module load samtools

fasta=$1
gsize=$2
opt=$3

script=$VGP_PIPELINE/stats/

#asm=${fasta/.fasta/}
asm_name=${fasta/.fasta/}
asm_name=`basename $asm_name`
asm=$asm_name

if ! [ -e $fasta.fai ]; then
	samtools faidx $fasta
fi
awk '{print $1"\t"$2}' $fasta.fai > $fasta.len

echo
echo "=== Scaffolds ==="
java -jar -Xmx2g $script/lenCalcNGStats.jar $fasta.len $gsize 1 >> $asm.stats
#cat $asm.stats
echo "" >> $asm.stats


if [[ "$opt" == "s" ]]; then
        exit 0
fi

if [ ! -e $asm.contigs.len ]; then
	if [ ! -e $asm_name.gaps ]; then
		java -jar -Xmx4g $script/fastaGetGaps.jar $fasta $asm_name.gaps
	fi
	awk -F "\t" '{print $1"\t"$2"\t"$3}' $asm_name.gaps > $asm_name.gaps.bed
	awk '{print $1"\t0\t"$2}' $fasta.fai > $fasta.len.bed
	num_gaps=`wc -l $asm_name.gaps.bed | awk '{print $1}'`
	if [[ $num_gaps -gt 0 ]]; then
		echo "$num_gaps gaps found. Compute contig stats..."
		module load bedtools
		bedtools subtract -a $fasta.len.bed -b $asm_name.gaps.bed | awk '{print $1"\t"$NF-$(NF-1)}' > $asm_name.contigs.len
	fi
	rm $asm_name.gaps.bed
fi

if [[ $num_gaps -gt 0 ]]; then
	echo
	echo "=== Contigs ==="
	#if [ ! -e $asm.contigs.stats ]; then
		java -jar -Xmx2g $script/lenCalcNGStats.jar $asm_name.contigs.len $gsize 1 > $asm.contigs.stats
	#fi
	#cat $asm.contigs.stats
	echo "" >> $asm.contigs.stats

	echo
	echo "=== Gaps ==="
	#if [ ! -e $asm.gaps.stats ]; then
		java -jar -Xmx2g $script/lenCalcNGStats.jar $asm_name.gaps $gsize 3 > $asm.gaps.stats
	#fi
	#cat $asm.gaps.stats
	echo "" >> $asm.gaps.stats
	rm $asm_name.contigs.len
fi
rm $asm_name.gaps

if [[ "$opt" != "p" ]]; then
	exit 0
fi
echo

echo "=== Primary Stats ==="
grep "scaffold_" $fasta.len > $asm.p.len
java -jar -Xmx2g $script/lenCalcNGStats.jar $asm.p.len $gsize > $asm.p.stats

echo "Scaffolds"
cat $asm.p.stats
echo

grep "scaffold_" $asm.contigs.len > $asm.contigs.p.len
java -jar -Xmx2g $script/lenCalcNGStats.jar $asm.contigs.p.len $gsize 1 > $asm.contigs.p.stats
echo "Contigs"
cat $asm.contigs.p.stats
echo

grep "scaffold_" $asm.gaps > $asm.gaps.p
java -jar -Xmx2g $script/lenCalcNGStats.jar $asm.gaps.p $gsize 3 > $asm.gaps.p.stats
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
java -jar -Xmx2g $script/lenCalcNGStats.jar $asm.h.len $gsize > $asm.h.stats
cat $asm.h.stats
echo

#:<<'END'
if [ ! -e $asm.h.fasta ]; then
    echo "Extract alt set"
    cut -f1 $asm.h.len > $asm.h.list
    java -jar -Xmx2g $script/fastaExtractFromList.jar $asm.fasta $asm.h.list $asm.h.fasta
fi
#END
