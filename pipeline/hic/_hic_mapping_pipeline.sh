#!/bin/bash

########################################################################
# HIC MAPPING PIPELINE BASED ON ARIMA GENOMICS MAPPING PIPELINE 042817 #
########################################################################

# Below find the commands currently used to map HiC data for use in genotyping analysis. 

# Replace the variables at the top with the correct paths for the locations of files/programs on your system. 

# This bash script will map one paired end sample (read1 & read2 fastqs). Feel to modify and multiplex as you see fit to work with  your volume of samples and system. 

# This script has been modified to fit in a generalized pipelines

# Dependencies: bwa, samtools, picard <- dedup step only

if [[ -z $1 ]] ; then
        echo "Usage: ./_hic_mapping_pipeline.sh <ref.fasta> <prefix> <fastq.map> [job dependency]"
        echo "fastq.map: /path/to/R1.fastq.gz   /path/to/R2.fastq.gz tab-delimited, one per line"
        exit -1
fi

REF=$1 # ref fasta file
PREFIX=$2 # prefix for out files
FASTQ_MAP=$3 # fastq path file: /path/to/R1.fastq /path/to/R2.fastq
THREADS=$4 # num. threads to run
PARTITION=$5 # slurm partition
if [ ! -z $5 ]; then
    PARTITION="-p $PARTITION"
fi
SCRATCH=$6 # scratch directory, default to /lscratch
if [ -z $6 ]; then
    SCRATCH="/lscratch"
fi
BWA_OPTS=$7
if [ -z "$7" ]; then
    BWA_OPTS="-B8"
fi

awk '{printf $1"\n"$2"\n"}' $FASTQ_MAP > fastq.ls

BWA='bwa'
FILTER="$VGP_PIPELINE/hic/filter_five_end.pl"
COMBINER="$VGP_PIPELINE/hic/two_read_bam_combiner.pl"

RAW_DIR="$SCRATCH/raw"
FILT_DIR="$SCRATCH/filtered"
COMB_DIR="$SCRATCH/combined"

mkdir -p $RAW_DIR $FILT_DIR $COMB_DIR logs

echo "Using VGP assembly pipeline in: $VGP_PIPELINE"

echo "### Step 0: FASTA indexing"
echo "\
sbatch $PARTITION -c1 -o logs/index.%A.out $VGP_PIPELINE/hic/index.sh $REF | awk '{print $4}' > index.jid"
sbatch $PARTITION -c1 -o logs/index.%A.out $VGP_PIPELINE/hic/index.sh $REF | awk '{print $4}' > index.jid

echo "\
sbatch $PARTITION --array=1-$(wc -l < fastq.ls) --dependency=afterok:$(cat index.jid) --cpus-per-task=$THREADS -o logs/map.%A_%a.out $VGP_PIPELINE/hic/map.sh $BWA $THREADS $BWA_OPTS $REF $RAW_DIR $PREFIX | awk '{print $4}' > map.jid"
sbatch $PARTITION --array=1-$(wc -l < fastq.ls) --dependency=afterok:$(cat index.jid) --cpus-per-task=$THREADS -o logs/map.%A_%a.out $VGP_PIPELINE/hic/map.sh $BWA $THREADS $BWA_OPTS $REF $RAW_DIR $PREFIX | awk '{print $4}' > map.jid

echo "\
sbatch $PARTITION --array=1-$(wc -l < fastq.ls) --dependency=aftercorr:$(cat map.jid) --cpus-per-task=$THREADS -o logs/filter.%A_%a.out $VGP_PIPELINE/hic/filter.sh $RAW_DIR $FILTER $THREADS $FILT_DIR $PREFIX | awk '{print $4}' > filter.jid"
sbatch $PARTITION --array=1-$(wc -l < fastq.ls) --dependency=aftercorr:$(cat map.jid) --cpus-per-task=$THREADS -o logs/filter.%A_%a.out $VGP_PIPELINE/hic/filter.sh $RAW_DIR $FILTER $THREADS $FILT_DIR $PREFIX | awk '{print $4}' > filter.jid

echo "\
sbatch $PARTITION --dependency=afterok:$(cat filter.jid) --cpus-per-task=$THREADS -o logs/combine.%A.out $VGP_PIPELINE/hic/combine.sh $COMBINER $FILT_DIR $PREFIX $THREADS $COMB_DIR"
sbatch $PARTITION --dependency=afterok:$(cat filter.jid) --cpus-per-task=$THREADS -o logs/combine.%A.out $VGP_PIPELINE/hic/combine.sh $COMBINER $FILT_DIR $PREFIX $THREADS $COMB_DIR

echo "#### All jobs submitted!"

:<<'END'
echo "### Start to dedup"

echo "\
dedup.sh $COMB_DIR/$PREFIX.bam $THREADS"
dedup.sh $COMB_DIR/$PREFIX.bam $THREADS
END


