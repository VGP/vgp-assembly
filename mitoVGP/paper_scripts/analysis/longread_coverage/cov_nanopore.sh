#!/bin/bash

HEADER=$(head -1 $1)
SEQ=$(tail -1 $1)
LN=$(tail -1 $1 | tr -d "\n" | wc -c)

printf "sequence length ${LN}"

printf "%s\n%s%s\n" $HEADER $SEQ $SEQ > dup_$1

seqtk subseq $2 reads.ls > filtered.fq

minimap2 -x map-ont -t 32 dup_$1 filtered.fq -a --secondary=no | samtools view -S -b -F 4 -F 0x800 | samtools sort > $3.bam

samtools depth -a $3.bam | awk '{print $3}' | head -${LN} > $3.head

samtools depth -a $3.bam | awk '{print $3}' | tail -${LN} > $3.tail

paste $3.head $3.tail | awk '{print $1+$2}' > $3.cov
