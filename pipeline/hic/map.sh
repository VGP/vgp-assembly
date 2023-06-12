#!/bin/bash

N=${SLURM_ARRAY_TASK_ID}
FASTQ=$(sed -n "${N}p" fastq.ls)

BWA=$1
THREADS=$2
BWA_OPTS=$3
REF=$4
RAW_DIR=$5
PREFIX=$6

echo "### Step 1.$N: FASTQ to BAM"
echo "\
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ | samtools view -Sb - > $RAW_DIR/${PREFIX}_${N}.bam"
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ | samtools view -Sb - > $RAW_DIR/${PREFIX}_${N}.bam