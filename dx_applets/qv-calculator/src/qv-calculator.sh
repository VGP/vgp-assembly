#!/bin/bash
# qv-calculator 0.0.1

set -x -e -o pipefail

main() {

	sudo chmod 777 /usr/bin/samtools
	sudo chmod 777 /usr/bin/bcftools

    echo "Value of aligned_bam: '$BAM'"
    echo "Value of variant: '$VAR'"

    bam_name=$BAM_name
    var_name=$VAR_name

    dx download "$BAM" -o ${bam_name}

    dx download "$VAR" -o ${var_name}

	if [ -e aligned.genomecov ]; then
		echo "aligned.genomecov already exists. skip..."
	else
		echo "Collect coverage"
		echo "\
		samtools view -F 0x100 -hb $bam | bedtools genomecov -ibam - -split > aligned.genomecov"
		samtools view -F 0x100 -hb ${bam_name} | bedtools genomecov -ibam - -split > aligned.genomecov
	fi


	echo "\
	awk '{if ($1=="genome" && $2>3) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp3"
	awk '{if ($1=="genome" && $2>3) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp3
	NUM_BP3=`cat $genome.numbp3`
	echo "Total bases > 3x: echo $NUM_BP3" >> qv_report.txt

	awk '{if ($1=="genome" && $2>5) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp5"
	awk '{if ($1=="genome" && $2>5) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp5
	NUM_BP5=`cat $genome.numbp5`
	echo "Total bases > 5x: echo $NUM_BP5" >> qv_report.txt

	bcftools view -H ${var_name} | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > $genome.numvar

	if [ -e $genome.numvar ]; then
		NUM_VAR=`cat $genome.numvar`
		echo "Total num. bases subject to change: $NUM_VAR" >> qv_report.txt
		QV3=`echo "$NUM_VAR $NUM_BP3" | awk '{print (-10*log($1/$2)/log(10))}'`
		echo $QV3 > $genome.qv3
		echo "Approximate QV of this genome $genome: $QV3" >> qv_report.txt

		QV5=`echo "$NUM_VAR $NUM_BP5" | awk '{print (-10*log($1/$2)/log(10))}'`
		echo $QV5 > $genome.qv5
		echo "QV of regiones above >5x for $genome: $QV5" >> qv_report.txt		
		
	fi


    qv_report=$(dx upload qv_report.txt --brief)

    dx-jobutil-add-output qv_report "$qv_report" --class=file
    
    genomecov=$(dx upload aligned.genomecov --brief)

    dx-jobutil-add-output genomecov "$genomecov" --class=file
}
