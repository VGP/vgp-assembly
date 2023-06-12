#!/bin/bash

N=${SLURM_ARRAY_TASK_ID}
FASTQ=$(sed -n "${N}p" fastq.ls)

RAW_DIR=$1
FILTER=$2
THREADS=$3
FILT_DIR=$4
PREFIX=$5

echo "### Step 2.$i: Filter 5' end"
echo "\
samtools view -h $RAW_DIR/${PREFIX}_${N}.bam | perl $FILTER | samtools view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_${N}.bam"
samtools view -h $RAW_DIR/${PREFIX}_${N}.bam | perl $FILTER | samtools view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_${N}.bam

echo "### Remove intermediate file"
echo "\
rm $RAW_DIR/${PREFIX}_${N}.bam"
#rm $RAW_DIR/${PREFIX}_${N}.bam