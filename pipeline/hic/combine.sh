#!/bin/bash

COMBINER=$1
FILT_DIR=$2
PREFIX=$3
THREADS=$4
COMB_DIR=$5

echo "### Step 3: Combine multiple input files"

N=$(ls $FILT_DIR/ | wc -l)
FW=$(printf "$FILT_DIR/${PREFIX}_%s.bam " $(seq 1 2 $N))
RV=$(printf "$FILT_DIR/${PREFIX}_%s.bam " $(seq 2 2 $N))

echo "\
samtools cat $FW -o $FILT_DIR/${PREFIX}_fw.bam"
samtools cat $FW -o $FILT_DIR/${PREFIX}_fw.bam

echo "\
samtools cat $RV -o $FILT_DIR/${PREFIX}_rv.bam"
samtools cat $RV -o $FILT_DIR/${PREFIX}_rv.bam

echo "### Step 4: Filter Combiner"
echo "\
perl $COMBINER $FILT_DIR/${PREFIX}_fw.bam $FILT_DIR/${PREFIX}_rv.bam | samtools view -@$THREADS -Sb > $COMB_DIR/$PREFIX.bam"
perl $COMBINER $FILT_DIR/${PREFIX}_fw.bam $FILT_DIR/${PREFIX}_rv.bam | samtools view -@$THREADS -Sb > $COMB_DIR/$PREFIX.bam

echo "rm $FILT_DIR/${PREFIX}*bam"
#rm $FILT_DIR/${PREFIX}*bam

echo "mv $COMB_DIR/$PREFIX.bam ."
mv $COMB_DIR/$PREFIX.bam .

echo "#### Finished Mapping!"