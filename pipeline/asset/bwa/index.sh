#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./index.sh <asm.fasta>"
	exit -1
fi

module load bwa/0.7.17

echo "
bwa index $1"
bwa index $1
