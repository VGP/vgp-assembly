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
	threads=4
fi

if [ ! -s $sample.bcf ] ; then
	ls --color=never bcf/*.bcf > concat_list.txt
	echo "\
	bcftools concat -f concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o $sample.bcf --threads $threads" &&
	bcftools concat -f concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o $sample.bcf --threads $threads &&
	echo || exit -1
fi

if [ ! -s $sample.bcf.csi ]; then
	echo "\
	bcftools index $sample.bcf" &&
	bcftools index $sample.bcf &&
	echo || exit -1
fi

echo "\
bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa") && INFO/DP>5 && (FORMAT/AD[:1]) / (FORMAT/AD[:1]+FORMAT/AD[:0]) > 0.5' -Oz --threads=$threads $sample.bcf > $sample.changes.vcf.gz" &&
bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa") && INFO/DP>5 && (FORMAT/AD[:1]) / (FORMAT/AD[:1]+FORMAT/AD[:0]) > 0.5' -Oz --threads=$threads $sample.bcf > $sample.changes.vcf.gz &&

echo "\
bcftools index $sample.changes.vcf.gz" &&
bcftools index $sample.changes.vcf.gz &&
echo &&

echo "\
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa") && INFO/DP>5' -Hla -f $fasta $sample.changes.vcf.gz > ${sample}_fb.fasta" &&
bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa") && INFO/DP>5' -Hla -f $fasta $sample.changes.vcf.gz > ${sample}_fb.fasta &&
echo &&

echo "Done!" ||
echo "Failed..."

echo
bcftools view -H -i 'QUAL>1 && (GT="AA" || GT="Aa") && INFO/DP>5' -Ov $sample.changes.vcf.gz | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > $sample.numvar
echo "Num. bases affected: `cat $sample.numvar`"
echo

