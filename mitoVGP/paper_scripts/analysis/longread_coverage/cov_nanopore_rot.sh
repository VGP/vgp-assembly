#!/bin/bash

HEADER=$(head -1 $1)
SEQ=$(tail -1 $1)
LN=$(tail -1 $1 | tr -d "\n" | wc -c)

CENTER=$(echo "scale = 0; (${4} + ${5}) / 2" | bc -l)
DIFF=$(echo "scale = 0; ${CENTER} - (${LN} / 2)" | bc -l)

echo $CENTER
echo $DIFF

END=$(tail -1 $1 | tr -d "\n" | cut -c1-${DIFF})
START=$(tail -1 $1 | tr -d "\n" | cut -c$((${DIFF} + 1))-${LN})

printf "%s\n%s%s\n" $HEADER $START $END > rot_$1

seqtk subseq $2 reads.ls > filtered.fq

minimap2 -x map-ont -t 32 rot_$1 filtered.fq -a --secondary=no | samtools view -S -b -F 4 -F 0x800 | samtools sort > $3.bam

samtools depth -a $3.bam | awk '{print $3}' | head -${LN} > $3.head

samtools depth -a $3.bam | awk '{print $3}' | tail -${LN} > $3.tail

paste $3.head $3.tail | awk '{print $1+$2}' > $3.cov
