#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./ast.sh <gaps.bed> [mean_cov]"
	echo "    When [mean_cov] is not provided, ast_10x will run twice to get the mean_cov"
	exit -1
fi

ls aligned.bam
gaps=$1
mean_cov=$2

platform=10x

if [[ -z $mean_cov ]]; then
	echo "First round to get mean_cov"

	if [[ -s $platform.bed ]]; then
		echo "*** 10x.bed found. skip the first round ***"
	else
		echo "
		$asset/bin/ast_10x -x $gaps aligned.bam > $platform.bed"
		$asset/bin/ast_10x -x $gaps aligned.bam > $platform.bed
		echo
	fi

	# Getting it from log: DEPRECATED as it is system dependent
	# grep coverage logs/10x.ast.*.log | tail -n1 | awk '{print $5}' | sed 's/,//g' | awk '{printf "%.0f\n", $1}'
	mean_cov=`awk '{sum+=$1*$2; total+=$2} END {printf "%.0f\n", sum/total}' TX.stat`
fi

echo "mean_cov : $mean_cov"
cutoff=`echo $mean_cov | awk '{printf "%.0f\n", $1*3.5}'`
echo "cutoff : $cutoff"

echo "Second round with high coverage threshold (-C) set"

if [[ ! -s ${platform}_C.bed ]]; then
	echo "
	$asset/bin/ast_10x -x -C $cutoff $gaps aligned.bam > ${platform}_C.bed"
	$asset/bin/ast_10x -x -C $cutoff $gaps aligned.bam > ${platform}_C.bed
	echo
fi

module load bedtools

bedtools subtract -a ../asm.bed -b ${platform}_C.bed | bedtools merge -d 100 -i - > $platform.low_high.bed
bedtools subtract -a $platform.low_high.bed -b ../asm.ends.bed -A > $platform.low_high.trim1k.bed
bedtools subtract -a ../asm.bed -b $platform.low_high.bed > $platform.support.bed

