#!/bin/bash

sample=$1
if [[ -z $sample ]]; then
    echo "sample=$sample"
    exit -1
fi

ref=${sample}_t1
fasta=refdata-$ref/fasta/genome.fa

module load samtools
#### [+] Loading samtools 1.8  ... 
# This is co-installed with bcftools

ls --color=never bcf/*.bcf > concat_list.txt
bcftools concat -nf concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o $sample.bcf
bcftools index $sample.bcf
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $fasta $sample.bcf > ${sample}_t2.fasta


