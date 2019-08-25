#!/bin/bash

module load bedtools
module load samtools
module load bcftools

if [ -z $1 ]; then
	echo "Usage: ./genomecov.sh <genome_id> [bam]"
	echo "Requires summary.csv, aligned.benomecov, and <genome_id>.changes.vcf.gz"
	exit -1
fi

bam=aligned.bam		#aligned, pos-sorted
genome=$1

if [ ! -z $2 ]; then
	bam=$2
fi

if [ -s aligned.genomecov ]; then
	echo "aligned.genomecov already exists. skip..."
else
	echo "Collect coverage"
	echo "\
	samtools view -F 0x100 -u $bam | bedtools genomecov -ibam - -split > aligned.genomecov"
	samtools view -F 0x100 -u $bam | bedtools genomecov -ibam - -split > aligned.genomecov
fi
echo


mean_cov=`tail -n1 summary.csv | awk -F "," '{printf "%.0f\n", $17}'`	# parse out the mean_cov from summary.csv
h_filter=$((mean_cov*12))	# exclude any sites >12x
l_filter=3			# exclude any sites <3x
echo "Get numbp between $l_filter ~ $h_filter x"

echo "\
awk -v l=$l_filter -v h=$h_filter '{if (\$1=="genome" && \$2>l && \$2<h) {numbp += \$3}} END {print numbp}' aligned.genomecov > $genome.numbp"
awk -v l=$l_filter -v h=$h_filter '{if ($1=="genome" && $2>l && $2<h) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp
NUM_BP=`cat $genome.numbp`
echo "Total bases > 3x: $NUM_BP"

if [ -e $genome.numvar ]; then
	bcftools view -H -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Ov $genome.changes.vcf.gz | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > $genome.numvar
	NUM_VAR=`cat $genome.numvar`
	echo "Total num. bases subject to change: $NUM_VAR"
	QV=`echo "$NUM_VAR $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}'`
	echo $QV > $genome.qv
	echo "QV of this genome $genome: $QV"
fi
