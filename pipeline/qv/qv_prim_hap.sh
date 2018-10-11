#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./qv_prim_hap.sh <genome_id> <primary_fasta> <alt_fasta>"
	exit -1
fi

genome=$1
fa_primary=$2
fa_alt=$3

threads=$SLURM_CPUS_PER_TASK
if [ -z $threads ]; then
	threads=2
fi

echo "Load modlues"
module load bedtools
module load samtools
echo

if [ ! -e $genome.changes.vcf.gz ]; then
	echo "Get changes"
	echo "\
	bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=$threads $genome.bcf > $genome.changes.vcf.gz"
	bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=$threads $genome.bcf > $genome.changes.vcf.gz
	echo "\
	bcftools index $genome.changes.vcf.gz"
	bcftools index $genome.changes.vcf.gz
	echo
fi

if [ ! -e $fa_primary.len ]; then
	echo "$fa_primary.len not found. Generateing..."
	java -jar -Xmx1g $VGP_PIPELINE/stats/fastaContigSize.jar $fa_primary
fi
awk '{print $1"\t0\t"$2}' $fa_primary.len > prim.bed

if [ ! -e $fa_alt.len ]; then
	echo "$fa_alt.len not found. Generating..."
	java -jar -Xmx1g $VGP_PIPELINE/stats/fastaContigSize.jar $fa_alt
fi
awk '{print $1"\t0\t"$2}' $fa_alt.len > alts.bed

for hap in prim alts
do
	echo "Collect $genome.$hap.numvar"
	bcftools view -H -R $hap.bed $genome.changes.vcf.gz | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > $genome.$hap.numvar
	echo "Num. bases affected: `cat $genome.$hap.numvar`"
	echo

	echo "Collect $hap contigs/scaffold names"
	cut -f1 $hap.bed > $hap.list
	echo

	echo "\
	java -jar -Xmx1g $VGP_PIPELINE/qv/txtContains.jar aligned.genomecov $hap.list 1 | awk '{if (\$2>3) {numbp+=\$3}} END {print numbp}' - > $genome.$hap.numbp"
	java -jar -Xmx1g $VGP_PIPELINE/qv/txtContains.jar aligned.genomecov $hap.list 1 | awk '{if ($2>3) {numbp+=$3}} END {print numbp}' - > $genome.$hap.numbp

	NUM_BP=`cat $genome.$hap.numbp`
	echo "Total bases > 3x: in $hap: $NUM_BP"
	NUM_VAR=`cat $genome.$hap.numvar`
	echo "Total num. bases subject to change in $hap: $NUM_VAR"
	QV=`echo "$NUM_VAR $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}'`
	echo $QV > $genome.$hap.qv
	echo "QV of this genome $genome $hap: $QV"
	echo
done
