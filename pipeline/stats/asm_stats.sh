#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./asm_stats.sh <asm.fasta> <exp_genome_size (bp)>"
    exit -1
fi

java -jar -Xmx1g /home/rhiea/codes/fastaContigSize.jar $1
java -jar -Xmx1g /home/rhiea/codes/lenCalcNGStats.jar $1.len 1031000000 > ${1/.fasta/.stats}
cat ${1/.fasta/.stats}
