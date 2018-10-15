#!/bin/bash

if [ -z $1 ] ; then
	echo "Assumes we have asm.fasta"
	exit -1
fi

ref=asm.fasta

echo "=== start indexing reference ==="
echo "\
$tools/longranger/longranger-2.2.2/longranger mkref $ref"
$tools/longranger/longranger-2.2.2/longranger mkref $ref
echo ""

