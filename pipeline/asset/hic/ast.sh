#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./ast.sh <gaps.bed>"
	exit -1
fi

ls -a *.bam

gaps=$1

if [[ ! -s hic.bed ]]; then
	echo "
	$asset/bin/ast_hic $gaps *.bam > hic.bed"
	$asset/bin/ast_hic $gaps *.bam > hic.bed
fi

module load bedtools

platform=hic
bedtools subtract -a ../asm.bed -b $platform.bed | bedtools merge -d 100 -i - > $platform.low_high.bed
bedtools subtract -a $platform.low_high.bed -b ../asm.ends.bed -A > $platform.low_high.trim1k.bed
bedtools subtract -a ../asm.bed -b $platform.low_high.bed > $platform.support.bed
