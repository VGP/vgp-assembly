#!/bin/sh

module load bwa/0.7.17
module load samtools/1.8
module load picard/2.9.2

fa=$1

echo "bwa index $fa"
bwa index $fa
echo "samtools faidx $fa"
samtools faidx $fa
echo "\
java -jar -Xmx2g /usr/local/apps/picard/2.9.2/build/libs/picard.jar CreateSequenceDictionary R=$fa O=${fa/.fasta/.dict}"
java -jar -Xmx2g /usr/local/apps/picard/2.9.2/build/libs/picard.jar CreateSequenceDictionary R=$fa O=${fa/.fasta/.dict}
