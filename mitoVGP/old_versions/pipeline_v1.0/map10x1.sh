#!/bin/bash

set -e

SPECIES=$1
ABBR=$2
CONTIG=$3

if ! [[ -e "${SPECIES}/${ABBR}_10x1" ]]; then

mkdir ${SPECIES}/${ABBR}_10x1

fi

if ! [[ -e "${SPECIES}/${ABBR}_10x1/${CONTIG}" ]]; then

cp ${SPECIES}/${ABBR}_arrow/${CONTIG} ${SPECIES}/${ABBR}_10x1/

fi

cd ${SPECIES}/${ABBR}_10x1

dw_date=`date "+%Y%m%d-%H%M%S"`; aws s3 --no-sign-request ls s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/10x/ > file_list_$dw_date.txt

if ! [[ -e "${ABBR}.1.bt2" ]]; then

bowtie2-build ${CONTIG} ${ABBR}
echo "--Index successfully built."

fi

if ! [[ -e "p_fw.txt" ]]; then

grep -o -E "${ABBR}.*_R1_.*" file_list_${dw_date}.txt | sort | uniq > p_fw.txt

fi

if ! [[ -e "p_rv.txt" ]]; then

grep -o -E "${ABBR}.*_R2_.*" file_list_${dw_date}.txt | sort | uniq > p_rv.txt

fi

mapfile -t p1 < p_fw.txt
mapfile -t p2 < p_rv.txt

echo "--Following PE files found:"

for ((i=0; i<${#p1[*]}; i++));
do

echo ${p1[i]} ${p2[i]} $i

done

if ! [[ -e "aligned_10x1_raw_reads" ]]; then

mkdir aligned_10x1_raw_reads

fi

for ((i=0; i<${#p1[*]}; i++));
do

if ! [[ -e "aligned_10x1_raw_reads/aligned_${ABBR}_${i}.bam" ]]; then

echo "--Retrieve and align:"

echo "s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/10x/${p1[i]}"
echo "s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/10x/${p2[i]}"

aws s3 --no-sign-request cp s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/10x/${p1[i]} .
aws s3 --no-sign-request cp s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/10x/${p2[i]} .

bowtie2 -x ${ABBR} -1 ${p1[i]} -2 ${p2[i]} -p 12 | samtools view -bSF4 - > "aligned_10x1_raw_reads/aligned_${ABBR}_${i}.bam"
rm ${p1[i]} ${p2[i]}

fi

done

if ! [[ -e "aligned_${ABBR}_all.bam" ]]; then

samtools merge aligned_${ABBR}_all.bam aligned_10x1_raw_reads/*.bam

fi

READ_N=$(bedtools bamtobed -i aligned_${ABBR}_all.bam | cut -f 4 | wc -l)

if (("$READ_N" >= "1500000")); then

N=$(awk "BEGIN {printf \"%.2f\",1500000/${READ_N}}")
N_READ_N=$(awk "BEGIN {printf \"%.0f\",$N*${READ_N}}")

echo "The number of reads is above 1.5M ($READ_N). Subsampling by $N ($N_READ_N)"

samtools view -s $N -b aligned_${ABBR}_all.bam > aligned_${ABBR}_all_sub.bam
mv aligned_${ABBR}_all.bam aligned_${ABBR}_all_o.bam
mv aligned_${ABBR}_all_sub.bam aligned_${ABBR}_all.bam

fi

if ! [[ -e "fq" ]]; then

samtools fastq aligned_${ABBR}_all.bam -1 aligned_${ABBR}_all_1.fq -2 aligned_${ABBR}_all_2.fq -s aligned_${ABBR}_all_s.fq
mkdir fq
mv *.fq fq/

fi

if ! [[ -e "aligned_${ABBR}_all_sorted.bam" ]]; then

samtools sort aligned_${ABBR}_all.bam -o aligned_${ABBR}_all_sorted.bam -@ 12
samtools index aligned_${ABBR}_all_sorted.bam

fi

#simple-consensus-per-read-group aligned_${ABBR}_all_sorted.bam --max-depth 100000 --ignore-read-groups | sed -n '1~4s/^@/>/p;2~4p' > ${ABBR}.contig${CONTIG}_arrow2_10x1.fasta