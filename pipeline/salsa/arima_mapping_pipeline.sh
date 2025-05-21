#!/bin/bash

##########################################
# ARIMA GENOMICS MAPPING PIPELINE 042817 #
##########################################

#Below find the commands currently used to map HiC data for use in genotyping analysis. 

#Replace the variables at the top with the correct paths for the locations of files/programs on your system. 

#This bash script will map one paired end sample (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with  your volume of samples and system. 

# This script has been modified to fit in a generalized pipelines

##########################################
# Commands #
##########################################

fastq_map=$1  # fastq path file: /path/to/R1.fastq <tab> /path/to/R2.fastq
PREFIX=$2     # prefix for out files
REF=$3        # ref fasta file with indexes
SCRATCH=$4    # scratch directory, default to /lscratch
if [ -z $4 ]; then
    SCRATCH="/lscratch"
fi
BWA_OPTS=$5
if [ -z "$5" ]; then
    BWA_OPTS="-B8"
fi

THREADS=$SLURM_CPUS_PER_TASK	# num. threads to run

FASTQ1=`awk '{print $1}' $fastq_map`
FASTQ2=`awk '{print $2}' $fastq_map`

module load bwa/0.7.17
module load samtools/1.21

BWA='bwa'
SAMTOOLS='samtools'
FILTER="$VGP_PIPELINE/salsa/filter_five_end.pl"
COMBINER="$VGP_PIPELINE/salsa/two_read_bam_combiner.pl"

RAW_DIR="$SCRATCH/$SLURM_JOBID/raw"
FILT_DIR="$SCRATCH/$SLURM_JOBID/filtered"
COMB_DIR="$SCRATCH/$SLURM_JOBID/combined"

mkdir -p $RAW_DIR
mkdir -p $FILT_DIR
mkdir -p $COMB_DIR

echo "### Step 1.A: FASTQ to BAM (1st)"
echo "\
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ1 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_1.bam"
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ1 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_1.bam
echo ""

echo "### Step 2.A: Filter 5' end (1st)"
echo "\
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_1.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_1.bam"
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_1.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_1.bam &&
echo "" || exit -1

echo "### Remove intermediate file"
echo "\
rm $RAW_DIR/${PREFIX}_1.bam"
rm $RAW_DIR/${PREFIX}_1.bam
echo

echo "### Step 1.B: FASTQ to BAM (2nd)"
echo "\
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ2 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_2.bam"
$BWA mem -t$THREADS $BWA_OPTS $REF $FASTQ2 | $SAMTOOLS view -Sb - > $RAW_DIR/${PREFIX}_2.bam
echo ""

echo "### Step 2.B: Filter 5' end (2nd)"
echo "\
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_2.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_2.bam"
$SAMTOOLS view -h $RAW_DIR/${PREFIX}_2.bam | perl $FILTER | $SAMTOOLS view -@$THREADS -Sb - > $FILT_DIR/${PREFIX}_2.bam &&
echo "" || exit -1

echo "### Remove intermediate file"
echo "\
rm $RAW_DIR/${PREFIX}_2.bam"
rm $RAW_DIR/${PREFIX}_2.bam
echo

echo "### Step 3.A: Filter Combiner"
echo "\
perl $COMBINER $FILT_DIR/${PREFIX}_1.bam $FILT_DIR/${PREFIX}_2.bam | $SAMTOOLS view -@$THREADS -Sb > $COMB_DIR/$PREFIX.bam"
perl $COMBINER $FILT_DIR/${PREFIX}_1.bam $FILT_DIR/${PREFIX}_2.bam | $SAMTOOLS view -@$THREADS -Sb > $COMB_DIR/$PREFIX.bam
echo ""

echo "mv $COMB_DIR/$PREFIX.bam ."
mv $COMB_DIR/$PREFIX.bam .

echo "#### Finished Mapping!"
echo ""

:<<'END'
echo "### Start dedup (optional)"

echo "\
dedup.sh $COMB_DIR/$PREFIX.bam $THREADS"
dedup.sh $COMB_DIR/$PREFIX.bam $THREADS
END


