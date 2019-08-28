#!/bin/bash

set -e

#set variables species abbreviation
ABBR=$1
CONTIG=$2


FNAME="${CONTIG%%.*}"
FEXT="${CONTIG#*.}"
CONTIG_NAME=$(cat ${CONTIG} | awk '$0 ~ ">" {print substr($0,2)}')

if ! [[ -e "realigned_${ABBR}_all_sorted.bam" ]]; then

samtools sort -n aligned_${ABBR}_all_sorted.bam -o aligned_${ABBR}_all_paired.bam
samtools fastq aligned_${ABBR}_all_paired.bam -1 aligned_${ABBR}_all_1.fq -2 aligned_${ABBR}_all_2.fq -s aligned_${ABBR}_all_s.fq
bowtie2-build ${CONTIG} ${ABBR}
bowtie2 -x ${ABBR} -1 aligned_${ABBR}_all_1.fq -2 aligned_${ABBR}_all_2.fq -p 16 --no-mixed | samtools view -bSF4 - > "realigned_${ABBR}_all.bam"
samtools sort realigned_${ABBR}_all.bam -o realigned_${ABBR}_all_sorted.bam -@ 16
samtools index realigned_${ABBR}_all_sorted.bam

fi

if ! [[ -e "${FNAME}.delta" ]]; then

nucmer --maxmatch --nosimplify ${CONTIG} ${CONTIG} -f --delta "${FNAME}.delta" -b 500 -t 16 

fi

NUCMER_OUT=$(show-coords "${FNAME}.delta" -lrcTHo | grep "BEGIN" | head -1)
BEGIN1=$(echo $NUCMER_OUT | awk '{print $1}')
BEGIN2=$(echo $NUCMER_OUT | awk '{print $2}')
END1=$(echo $NUCMER_OUT | awk '{print $3}')

if (( ${BEGIN2} > ${END1} )); then

echo ">${FNAME}" > "${FNAME}_new.fasta" & grep -v ">" ${CONTIG} | tr -d '\n' | cut -c${BEGIN1}-${BEGIN2} >> "${FNAME}_new.fasta"

nucmer --maxmatch --nosimplify "${FNAME}_new.fasta" "${FNAME}_new.fasta" -f --delta "${FNAME}.delta" -b 500 -t 16 

CONTIG="${FNAME}_new.fasta"

NUCMER_OUT=$(show-coords "${FNAME}.delta" -lrcTHo | grep "BEGIN" | head -1)
BEGIN1=$(echo $NUCMER_OUT | awk '{print $1}')
BEGIN2=$(echo $NUCMER_OUT | awk '{print $2}')
END1=$(echo $NUCMER_OUT | awk '{print $3}')
	
fi

END2=$(echo $NUCMER_OUT | awk '{print $4}')
MIDDLE="$(( ${BEGIN2} + 1 ))-$(( ${END1} - 1 ))"

echo ">${FNAME}" > "${FNAME}_final.fasta" & grep -v ">" ${CONTIG} | tr -d '\n' | cut -c${MIDDLE} >> "${FNAME}_final.fasta"

echo ">${FNAME}_begin_${BEGIN1}-${BEGIN2}" > "${FNAME}_ends.fasta" & grep -v ">" ${CONTIG} | tr -d '\n' | cut -c${BEGIN1}-${BEGIN2} >> "${FNAME}_ends.fasta"
echo ">${FNAME}_end_${END1}-${END2}" >> "${FNAME}_ends.fasta" & grep -v ">" ${CONTIG} | tr -d '\n' | cut -c${END1}-${END2} >> "${FNAME}_ends.fasta"

arrCOV1=($(samtools depth -aa -r ${CONTIG_NAME}:${BEGIN1}-${BEGIN2} --reference ${CONTIG} realigned_${ABBR}_all_sorted.bam | awk '{print $3}'))
arrCOV2=($(samtools depth -aa -r ${CONTIG_NAME}:${END1}-${END2} --reference ${CONTIG} realigned_${ABBR}_all_sorted.bam | awk '{print $3}'))

if ! [[ -e "${FNAME}_ends_aligned.fasta" ]]; then

muscle -in ${FNAME}_ends.fasta -out ${FNAME}_ends_aligned.fasta

fi

if ! [[ -e "${FNAME}_ends_aligned.table" ]]; then

NAM1=$(cat ${FNAME}_ends_aligned.fasta | awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' | awk 'NR==1')
SEQ1=$(cat ${FNAME}_ends_aligned.fasta | awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' | awk 'NR==2')
NAM2=$(cat ${FNAME}_ends_aligned.fasta | awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' | awk 'NR==3')
SEQ2=$(cat ${FNAME}_ends_aligned.fasta | awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' | awk 'NR==4')

arrSEQ1=($(fold -w1 <<< "$SEQ1"))
arrSEQ2=($(fold -w1 <<< "$SEQ2"))

printf "%s\n" ${NAM1} > ${FNAME}_ends_aligned_1.table
printf "%s\n" "${arrSEQ1[@]}" >> ${FNAME}_ends_aligned_1.table

COUNTER=0

while read i;
	do
		if [ "$i" == ${NAM1} ]; then
			printf "%s\t%s\n" ${NAM1} "Cov" > ${FNAME}_ends_aligned_1_with_cov.table
		elif [ "$i" == "-" ]; then
			printf "%s\t%s\n" "$i" "NA" >> ${FNAME}_ends_aligned_1_with_cov.table
		else
			printf "%s %s\n" "$i" "${arrCOV1[$COUNTER]}" >> ${FNAME}_ends_aligned_1_with_cov.table
			let COUNTER=COUNTER+1	
		fi
	done < ${FNAME}_ends_aligned_1.table

printf "%s\n" $NAM2 > ${FNAME}_ends_aligned_2.table
printf "%s\n" "${arrSEQ2[@]}" >> ${FNAME}_ends_aligned_2.table

COUNTER=0

while read i;
	do
		if [ "$i" == ${NAM2} ]; then
			printf "%s\t%s\n" ${NAM2} "Cov" > ${FNAME}_ends_aligned_2_with_cov.table
		elif [ "$i" == "-" ]; then
			printf "%s\t%s\n" "$i" "NA" >> ${FNAME}_ends_aligned_2_with_cov.table
		else
			printf "%s %s\n" "$i" "${arrCOV2[$COUNTER]}" >> ${FNAME}_ends_aligned_2_with_cov.table
			let COUNTER=COUNTER+1
		fi
	done < ${FNAME}_ends_aligned_2.table

paste ${FNAME}_ends_aligned_1_with_cov.table ${FNAME}_ends_aligned_2_with_cov.table > ${FNAME}_ends_aligned.table

fi

S=""

COUNTER=0

arrN1=( $(awk 'FNR == 1 {next} {print $1}' ${FNAME}_ends_aligned.table) )
arrQ1=( $(awk 'FNR == 1 {next} {print $2}' ${FNAME}_ends_aligned.table) )
arrN2=( $(awk 'FNR == 1 {next} {print $3}' ${FNAME}_ends_aligned.table) )
arrQ2=( $(awk 'FNR == 1 {next} {print $4}' ${FNAME}_ends_aligned.table) )

while [  $COUNTER -lt ${#arrN1[@]} ]
	do
		N1=${arrN1[$COUNTER]}
		Q1=${arrQ1[$COUNTER]}
		N2=${arrN2[$COUNTER]}
		Q2=${arrQ2[$COUNTER]}
		VER=$COUNTER
		if [  ${N1} == ${N2} ]; then
			S="${S}${N1}"
		elif [  ${N1} == "-" ]; then
			while [  ${arrN1[${VER}]} == "-" ]; do
				if [[ ${arrQ1[$VER-1]} == ?(-)+([0-9]) ]] && (( ${arrQ1[$VER-1]} > ${Q2} )); then
				:
				elif [[ ${arrQ1[$VER-1]} == ?(-)+([0-9]) ]] && (( ${arrQ1[$VER-1]} < ${Q2} )); then
					S="${S}${N2}"
				fi			
    			let VER-=1
         	done
		elif [  ${N2} == "-" ]; then
			while [  ${arrN2[$VER]} == "-" ]; do
				if [[ ${arrQ2[$VER-1]} == ?(-)+([0-9]) ]] && (( ${arrQ2[$VER-1]} > ${Q1} )); then
				:
				elif [[ ${arrQ2[$VER-1]} == ?(-)+([0-9]) ]] && (( ${arrQ2[$VER-1]} < ${Q1} )); then
					S="${S}${N1}"
				fi			
    			let VER-=1
         	done
		elif [  ${N1} != "-" ] && [  ${N1} != ${N2} ]; then
				if ((  ${Q1} > ${Q2} )); then
					S="${S}${N1}"
				else
					S="${S}${N2}"
				fi				
		fi
		
		let COUNTER=COUNTER+1
		
	done

sed -i "$ s/$/$S/" ${FNAME}_final.fasta