#!/bin/bash

HEADER=$(head -1 $1)
SEQ=$(tail -1 $1)
LN=$(tail -1 $1 | tr -d "\n" | wc -c)

pbmm2 align -j 24 --best-n 1 $1 filtered.fq --sort | samtools view -S -b > $2.bam

samtools index $2.bam

samtools depth -a $2.bam | awk '{print $3}' | head -${LN} > $2.head

samtools depth -a $2.bam | awk '{print $3}' | tail -${LN} > $2.tail

paste $2.head $2.tail | awk '{print $1+$2}' > $2.cov
