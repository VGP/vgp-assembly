#!/bin/bash

HEADER=$(head -1 $1)
SEQ=$(tail -1 $1)
LN=$(tail -1 $1 | tr -d "\n" | wc -c)

printf "%s\n%s%s" $HEADER $SEQ $SEQ > dup_$1

picard RevertSam I=$2 O=unmapped_$2 MAX_DISCARD_FRACTION=0.000 ATTRIBUTE_TO_CLEAR=XT ATTRIBUTE_TO_CLEAR=XN ATTRIBUTE_TO_CLEAR=AS ATTRIBUTE_TO_CLEAR=OC ATTRIBUTE_TO_CLEAR=OP SORT_ORDER=unsorted RESTORE_ORIGINAL_QUALITIES=true REMOVE_DUPLICATE_INFORMATION=true REMOVE_ALIGNMENT_INFORMATION=true VALIDATION_STRINGENCY=LENIENT

pbmm2 align -j 24 --best-n 1 dup_$1 unmapped_$2 --sort | samtools view -S -b > realigned_$3.bam

samtools depth -a realigned_$3.bam | awk '{print $3}' | head -${LN} > $3.head

samtools depth -a realigned_$3.bam | awk '{print $3}' | tail -${LN} > $3.tail

paste $3.head $3.tail | awk '{print $1+$2}' > $3.cov
