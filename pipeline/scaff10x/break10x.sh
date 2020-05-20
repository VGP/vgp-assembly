#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./break10x.sh <out_prefix>"
    echo -e "\tSymlink asm.fasta, read-BC1.fastq and read-BC2.fastq"
    exit -1
fi

asm=asm.fasta

if [ ! -e asm.fasta ]; then
    echo "No asm.fasta detected. Symlink asm.fasta"
    exit -1
fi

if [ ! -e read-BC1.fastq ]; then
    echo "No read-BC1.fastq detected."
    exit -1
fi

out=$1  # Output prefix

$tools/scaff10x/Scaff10X_NoSam/Scaff10X/src/break10x \
    -nodes $SLURM_CPUS_PER_TASK -gap 100 \
    -reads 5 -score 20 \
    -cover 50 -ratio 15 \
    asm.fasta \
    read-BC1.fastq read-BC2.fastq \
    $out.fasta $out.breaks


