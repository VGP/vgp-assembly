#!/bin/bash

sample=$1
if [[ -z $sample ]]; then
    echo "sample=$sample"
    exit -1
fi

ref=asm
fasta=refdata-$ref/fasta/genome.fa

module load samtools
#### [+] Loading samtools 1.8  ... 
# This is co-installed with bcftools

threads=$SLURM_CPUS_PER_TASK
if [ -z $threads ]; then
	threads=2
fi

ls --color=never bcf/*.bcf > concat_list.txt
echo "\
bcftools concat -f concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o $sample.bcf --threads $threads" &&
bcftools concat -f concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o $sample.bcf --threads $threads &&
echo &&

echo "\
bcftools index $sample.bcf" &&
bcftools index $sample.bcf &&
echo &&

echo "\
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $fasta $sample.bcf > ${sample}_fb.fasta" &&
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $fasta $sample.bcf > ${sample}_fb.fasta &&
echo &&

echo "Done!" ||
echo "Failed..."

echo
bcftools view -H -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Ov $sample.bcf | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > $sample.numvar
echo "Num. bases affected: `cat $sample.numvar`"
