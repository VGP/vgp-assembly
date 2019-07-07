#!/bin/bash
# qv-calculator 0.0.1

set -x -e -o pipefail

main() {

	sudo chmod 777 /usr/bin/samtools
	sudo chmod 777 /usr/bin/bcftools

    echo "Value of aligned_bam: '$BAM'"
    echo "Value of variant: '$VAR'"
    echo "Value of csv: '$CSV'"

    bam_name=$BAM_name
    var_name=$VAR_name
    csv_name=$CSV_name

    dx download "$BAM" -o ${bam_name}

    dx download "$VAR" -o ${var_name}
    
    dx download "$CSV" -o ${csv_name}

	if [ -e aligned.genomecov ]; then
		echo "aligned.genomecov already exists. skip..."
	else
		echo "Collect coverage"
		echo "\
		samtools view -F 0x100 -hb $bam | bedtools genomecov -ibam - -split > aligned.genomecov"
		samtools view -F 0x100 -hb ${bam_name} | bedtools genomecov -ibam - -split > aligned.genomecov
	fi


	mean_cov=`tail -n1 ${csv_name} | awk -F "," '{printf "%.0f\n", $17}'`	# parse out the mean_cov from summary.csv
	h_filter=$((mean_cov*12))	# exclude any sites >12x
	l_filter=5			# exclude any sites <5x
	echo "Get numbp between $l_filter ~ $h_filter x"

	echo "\
	awk -v l=$l_filter -v h=$h_filter '{if (\$1=="genome" && \$2>l && \$2<h) {numbp += \$3}} END {print numbp}' aligned.genomecov > $genome.numbp"
	awk -v l=$l_filter -v h=$h_filter '{if ($1=="genome" && $2>l && $2<h) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp
	NUM_BP=`cat $genome.numbp`
	echo "Total bases > 5x: $NUM_BP" >> qv_report.txt

	bcftools view -H -i 'QUAL>1 && (GT="AA" || GT="Aa") && INFO/DP>5' -Ov ${var_name} | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > $genome.numvar

	if [ -e $genome.numvar ]; then
		NUM_VAR=`cat $genome.numvar`
		echo "Total num. bases subject to change: $NUM_VAR"  >> qv_report.txt
		QV=`echo "$NUM_VAR $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}'`
		echo $QV > $genome.qv
		echo "QV of this genome $genome: $QV"  >> qv_report.txt
	fi

    qv_report=$(dx upload qv_report.txt --brief)

    dx-jobutil-add-output qv_report "$qv_report" --class=file
    
    genomecov=$(dx upload aligned.genomecov --brief)

    dx-jobutil-add-output genomecov "$genomecov" --class=file
}
