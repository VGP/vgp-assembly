#!/bin/bash

module load python
#### Python 2.7.15 :: Anaconda custom (64-bit)
module load perl/5.18.2
#### Loading Perl 5.18.2  ... 
module load R
#### Loading gcc  7.2.0  ... 
#### Loading GSL 2.4 for GCC 7.2.0 ... 
#### Loading openmpi 3.0.0  for GCC 7.2.0 
#### Loading R 3.5.0_build2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

scripts=$tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/scripts

ENZYME1=$1
ENZYME2=$2

perl $scripts/fa2cmap_multi_color.pl -h

if [ -z $ENZYME1 ]; then
	echo "Usage: fa_to_cmap.sh <enzyme1> [enzyme2]"
	echo "Assumes we have asm.fasta"
	exit -1
fi

if [ -z $ENZYME2 ]; then
	echo "\
	perl $scripts/fa2cmap_multi_color.pl -i asm.fasta -e $ENZYME1 1 -o $ENZYME1"
	perl $scripts/fa2cmap_multi_color.pl -i asm.fasta -e $ENZYME1 1 -o $ENZYME1
else
	echo "\
	perl $scripts/fa2cmap_multi_color.pl -i asm.fasta -e $ENZYME1 1 -e $ENZYME2 2 -o ${ENZYME1}_${ENZYME2}"
	perl $scripts/fa2cmap_multi_color.pl -i asm.fasta -e $ENZYME1 1 -e $ENZYME2 2 -o ${ENZYME1}_${ENZYME2}
fi

