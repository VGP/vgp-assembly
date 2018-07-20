#!/bin/bash

BMAP1=DLE1.cmap # DLE1.cmap
BMAP2=BSSSI.cmap # BSSSI.cmap
ENZYME1=DLE1	# DLE1
ENZYME2=BSSSI	# BSSSI
ASM=asm.fasta	# ln -s to the asm.fasta

mkdir -p output

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
        /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/TGH/hybridScaffold_two_enzymes_DLE1.xml \
        --BNGPath1 $BMAP1 --BNGPath2 $BMAP2 \
        --NGSPath $ASM \
        --OutputDir $PWD/output \
        --RefAlignerPath /data/Phillippy/tools/bionano/Solve3.2.1_04122018/RefAligner/7437.7523rel/avx/RefAligner \
        --Enzyme1 $ENZYME1 \
        --Enzyme2 $ENZYME2 \
	> scaffold.out
"

Rscript /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/runTGH.R \
        /data/Phillippy/tools/bionano/Solve3.2.1_04122018/HybridScaffold/04122018/TGH/hybridScaffold_two_enzymes_DLE1.xml \
        --BNGPath1 $BMAP1 --BNGPath2 $BMAP2 \
        --NGSPath $ASM \
        --OutputDir $PWD/output \
        --RefAlignerPath /data/Phillippy/tools/bionano/Solve3.2.1_04122018/RefAligner/7437.7523rel/avx/RefAligner \
        --Enzyme1 $ENZYME1 \
        --Enzyme2 $ENZYME2 \
        > scaffold.out 2>&1
