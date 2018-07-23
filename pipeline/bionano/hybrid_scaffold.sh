#!/bin/bash

name=$1
ENZYME=$2
BMAP=$ENZYME.cmap # DLE1.cmap
ASM=asm.fasta	# ln -s to the asm.fasta
CONFIG=$3
RefAligner=/data/Phillippy/tools/bionano/Solve3.2.1_04122018/RefAligner/7437.7523rel/avx/RefAligner

module load python
#### Python 2.7.15 :: Anaconda custom (64-bit)
module load perl/5.18.2
#### Loading Perl 5.18.4  ... 
module load R
#### Loading gcc  7.2.0  ... 
#### Loading GSL 2.4 for GCC 7.2.0 ... 
#### Loading openmpi 3.0.0  for GCC 7.2.0 
#### Loading R 3.5.0_build2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

echo "\
perl /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/hybridScaffold.pl \
        -n $ASM \
	-b $BMAP \
	-c $CONFIG \
	-r $RefAligner \
	-B 2 \
	-N 2 \
        -o $PWD/$name
"
perl /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/hybridScaffold.pl \
        -n $ASM \
	-b $BMAP \
	-c $CONFIG \
	-r $RefAligner \
	-B 2 \
	-N 2 \
        -o $PWD/$name
