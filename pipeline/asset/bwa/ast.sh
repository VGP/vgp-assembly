#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./ast.sh <gaps.bed>"
	exit -1
fi

ls *.bam

gaps=$1

echo "
$asset/bin/ast_10x $gaps *.bam > 10x.bed"
$asset/bin/ast_10x $gaps *.bam > 10x.bed
