#!/bin/bash

name=$1
BMAP1=$2.cmap # BSPQI.cmap
BMAP2=$3.cmap # BSSSI.cmap
ENZYME1=$2	# BSPQI
ENZYME2=$3	# BSSSI
ASM=asm.fasta	# ln -s to the asm.fasta
CONFIG=$4
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

echo "\
Rscript /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/runTGH.R \
        $CONFIG \
        --BNGPath1 $BMAP1 --BNGPath2 $BMAP2 \
        --NGSPath $ASM \
        --OutputDir $PWD/$name \
        --RefAlignerPath $RefAligner \
        --Enzyme1 $ENZYME1 \
        --Enzyme2 $ENZYME2
"

Rscript /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/runTGH.R \
        $CONFIG \
        --BNGPath1 $BMAP1 --BNGPath2 $BMAP2 \
        --NGSPath $ASM \
        --OutputDir $PWD/$name \
        --RefAlignerPath $RefAligner \
        --Enzyme1 $ENZYME1 \
        --Enzyme2 $ENZYME2
