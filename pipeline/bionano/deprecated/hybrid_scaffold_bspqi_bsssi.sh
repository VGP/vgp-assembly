#!/bin/bash

BMAP1=BSPQI.cmap # DSPQI.cmap
BMAP2=BSSSI.cmap # BSSSI.cmap
ENZYME1=BSPQI	# BSPQI
ENZYME2=BSSSI	# BSSSI
ASM=asm.fasta	# ln -s to the asm.fasta
CONFIG=/data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/TGH/hybridScaffold_two_enzymes.xml
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
Rscript /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/runTGH.R \
        $CONFIG \
        --BNGPath1 $BMAP1 --BNGPath2 $BMAP2 \
        --NGSPath $ASM \
        --OutputDir $PWD/$name \
        --RefAlignerPath $RefAligner \
        --Enzyme1 $ENZYME1 \
        --Enzyme2 $ENZYME2 \
	> scaffold.out
"

Rscript /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/runTGH.R \
        $CONFIG \
        --BNGPath1 $BMAP1 --BNGPath2 $BMAP2 \
        --NGSPath $ASM \
        --OutputDir $PWD/$name \
        --RefAlignerPath $RefAligner \
        --Enzyme1 $ENZYME1 \
        --Enzyme2 $ENZYME2 \
        > scaffold.out 2>&1
