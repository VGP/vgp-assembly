#!/bin/bash

assembly=asm.fasta
R1=read-BC1.fastq
R2=read-BC2.fastq
name=$1

/data/Phillippy/tools/scaff10x/Scaff10X/src/scaff10x

module load bwa

mkdir -p round1
mkdir -p round2

# Round 1
cd round1
if [ -e "round1.fasta" ]; then
	echo "round1.fasta already exists. Skip the first step..."
else
	ln -s ../$assembly
	ln -s ../$R1
	ln -s ../$R2

	/data/Phillippy/tools/scaff10x/Scaff10X/src/scaff10x -nodes $SLURM_CPUS_PER_TASK -longread 1 -gap 100 -matrix 2000 -reads 12 -link 10 -block 50000 $assembly $R1 $R2 round1.fasta
fi
assembly=round1.fasta

# Round 2
cd ../round2
ln -s ../round1/$assembly
ln -s ../$R1
ln -s ../$R2
/data/Phillippy/tools/scaff10x/Scaff10X/src/scaff10x -nodes $SLURM_CPUS_PER_TASK -longread 1 -gap 100 -matrix 2000 -reads 8 -link 10 -block 50000 $assembly $R1 $R2 round2.fasta

cd ..
ln -s round2/round2.fasta ${name}_s1.fasta
