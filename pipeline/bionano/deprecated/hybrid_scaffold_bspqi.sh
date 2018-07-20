#!/bin/bash

BMAP=DLE1.cmap # DLE1.cmap
ASM=asm.fasta	# ln -s to the asm.fasta
CONFIG=/data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/hybridScaffold_config.xml
RefAligner=/data/Phillippy/tools/bionano/Solve3.2.1_04122018/RefAligner/7437.7523rel/avx/RefAligner
name=$1

module load python
#### Python 2.7.15 :: Anaconda custom (64-bit)
module load perl/5.18.2
#### Loading Perl 5.18.4  ... 
module load R
#### Loading gcc  7.2.0  ... 
#### Loading GSL 2.4 for GCC 7.2.0 ... 
#### Loading openmpi 3.0.0  for GCC 7.2.0 
#### Loading R 3.5.0_build2

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
        -o $PWD/$name \
        > scaffold.out 2>&1
