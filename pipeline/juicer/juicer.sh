#!/bin/bash

module load juicer/1.5.6
module load samtools
module load bwa

if [ ! -e reference/chr.sizes ]; then
   cd reference
   bwa index asm.fasta
   python /usr/local/apps/juicer/juicer-1.5.6/misc/generate_site_positions.py MboI asm asm.fasta
   samtools faidx asm.fasta
   awk '{print $1"\t"$2}' asm.fasta.fai >  chr.sizes
   cd -
fi

/usr/local/apps/juicer/juicer-1.5.6/scripts/juicer.sh -z  `pwd`/reference/asm.fasta -y  `pwd`/reference/asm_MboI.txt  -D /usr/local/apps/juicer/juicer-1.5.6/ -d `pwd` -p `pwd`/reference/chr.sizes
