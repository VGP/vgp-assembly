#!/bin/bash

#rm -fr ${2}"_"${4}

start=${8}
end=${9}
type=${4}

mkdir -p ${2}"_"${type}

HEADER=$(head -1 $1)
SEQ=$(tail -1 $1)
LN=$(tail -1 $1 | tr -d "\n" | wc -c)

CENTER=$(echo "scale = 0; (${start} + ${end}) / 2" | bc -l)
DIFF=$(echo "scale = 0; ${CENTER} - (${LN} / 2)" | bc -l)

if (( $DIFF > 0 )); then

	END=$(tail -1 $1 | tr -d "\n" | cut -c1-${DIFF})
	START=$(tail -1 $1 | tr -d "\n" | cut -c$((${DIFF} + 1))-${LN})

else

	END=$(tail -1 $1 | tr -d "\n" | cut -c1-$((${LN} - ${DIFF#-})))
	START=$(tail -1 $1 | tr -d "\n" | cut -c$((${LN} - ${DIFF#-} + 1))-${LN})

fi

cd ${2}"_"${type}

printf "%s\n%s%s\n" $HEADER $START $END > rot_$(basename $1)

dw_date=`date "+%Y%m%d-%H%M%S"`;

rm -f long_reads_file_list_*

printf "Collecting data using: aws s3 ls s3://genomeark/species/${3}/${2}/genomic_data/pacbio/\n\n"
aws s3 --no-sign-request ls s3://genomeark/species/${3}/${2}/genomic_data/pacbio/ | grep -oP "m.*.subreads.bam" | uniq > long_reads_file_list_$dw_date.txt

if ! [[ -e "*subreads.bam" ]] && ! [[ -e "reads.fastq" ]]; then

	aws s3 --no-sign-request cp --recursive --exclude="*scrap*" --exclude="*xml*" --exclude="*bai*" --exclude="*pbi*" --include="*.subreads.bam" s3://genomeark/species/${3}/${2}/genomic_data/pacbio/ .

fi

printf "\n--Following long read files found:\n"

while read p; do

echo ${p}

done < long_reads_file_list_$dw_date.txt

printf "\n"

if ! [[ -e "reads.fastq" ]]; then

#for each Pacbio raw data file do
while read p; do

	if ! [[ -e "aligned_$(basename -- "${p%.*}").bam" ]] && ! [[ $p == *scraps* ]] && ! [[ $p == *.pbi ]] && ([[ $p == *.bam ]] || [[ $p == *.fastq ]] || [[ $p == *.fasta ]] || [[ $p == *.fa ]] || [[ $p == *.fq ]]); then

		pbmm2 align --best-n 1 --min-concordance-perc 0 rot_$(basename $1) $p aligned_${p%.*}.bam -j 24

		rm -f ${p}

	fi
	
done < long_reads_file_list_$dw_date.txt

for f in aligned_*.bam; do
	
	filename=$(basename -- "$f")
	filename="${filename%.*}"
	
	if ! [[ -e "${f}.pbi" ]]; then

		pbindex ${f}

	fi
	
	printf "convert: ${f} to fastq\n"
		
	bam2fastq ${f} -o "${filename}"
	mv ${filename} ${filename}.int.fastq
	rm $f
	rm ${f}.pbi
	
done

fi

if ! [[ -e "reads.fastq" ]]; then

	cat *.int.fastq > reads.fastq
	rm -f *.int.fastq

fi

seqtk seq -a reads.fastq > reads.fasta

makeblastdb -in reads.fasta -parse_seqids -dbtype nucl -out blast.db

blastn -dust no -outfmt "6 sseqid slen qcovs" -query ../$1 -db blast.db -max_target_seqs 100000 | sort -k2 -nr | uniq |awk -v gsize="${LN}" '{printf $0 "\t" gsize/$2*$3 "\n"}' | awk -v len1="${6}" -v len2="${7}" -v per1="${5}" -v per2="${10}" '{if($4>per1 && $2<len1 && $2>len2 && $3>per2) printf $0 "\n"}' > ${2}_reads.out

awk '{if ($4>70) print $1}' ${2}_reads.out > ${2}_reads.ls

seqtk subseq *fastq* ${2}_reads.ls > filtered.fq

pbmm2 align --best-n 1 rot_$(basename $1) filtered.fq --sort | samtools view -S -b > $2.bam

samtools index $2.bam

new_start=$(echo "scale = 0; (${LN} / 2) - ( (${end} - ${start}) / 2)" | bc -l)

new_end=$(echo "scale = 0; (${LN} / 2) + ( (${end} - ${start}) / 2)" | bc -l)

printf "%s\t%s\t%s\n" $(echo $HEADER | tr -d ">") $new_start $new_end > coords.bed

bedtools intersect -F 1.0 -b coords.bed -a $2.bam > region.bam

samtools sort region.bam -o sorted_region.bam

samtools index sorted_region.bam

java -jar /bin/jvarkit/dist/pcrclipreads.jar -B coords.bed sorted_region.bam |samtools view -q 1 -F 4 -Sbu - |samtools sort -o clipped.bam - && samtools index clipped.bam

java -jar /bin/jvarkit/dist/biostar84452.jar clipped.bam -o hard_clipped.bam

ref=$(awk '{print $3-$2}' coords.bed)

echo ${2}"_"${type} > ${2}"_"${type}_het.out

samtools view hard_clipped.bam | awk -v ref=$ref '{print length($10)-ref}' >> ${2}"_"${type}_het.out


