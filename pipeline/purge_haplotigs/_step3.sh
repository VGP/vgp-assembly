#!/bin/bash

module load purge_haplotigs/0~20180710.f4fd019

genome=$1

if ! [ -e asm.fasta ]; then
	ln -s ../../${genome}_c1.fasta asm.fasta
fi

if ! [ -e asm.fasta.faidx ]; then
	module load samtools
	samtools faidx asm.fasta
fi

threads=$SLURM_CPUS_PER_TASK

echo "STEP 3. Purge!"
echo "\
purge_haplotigs  purge  -g asm.fasta  -o ${genome}_curated  -c coverage_stats.csv  -b aligned.bam -t $threads -windowmasker"
purge_haplotigs  purge  -g asm.fasta  -o ${genome}_curated  -c coverage_stats.csv  -b aligned.bam -t $threads -windowmasker &&
echo "Done!" ||
echo "Failed..."
