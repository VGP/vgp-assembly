#!/bin/bash

SPECIES=$1
ID=$2
REF=$3

rm -r $ID

mkdir -p $ID
mkdir -p ${ID}/log/

dw_date=`date "+%Y%m%d-%H%M%S"`
	
printf "Collecting data using: aws s3 --no-sign-request ls s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/\n"
aws s3 --no-sign-request ls s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/ > ${ID}/log/short_reads_file_list_aws_$dw_date.txt
	
awk '{print $4}' ${ID}/log/short_reads_file_list_aws_$dw_date.txt > ${ID}/log/short_reads_file_list_$dw_date.txt

grep -o -E ".*R1.*" ${ID}/log/short_reads_file_list_${dw_date}.txt | sort | uniq > ${ID}/log/p_fw.txt
grep -o -E ".*R2.*" ${ID}/log/short_reads_file_list_${dw_date}.txt | sort | uniq > ${ID}/log/p_rv.txt

mapfile -t p1 < ${ID}/log/p_fw.txt
mapfile -t p2 < ${ID}/log/p_rv.txt

printf "\n--Following PE files found:\n"

for ((i=0; i<${#p1[*]}; i++));
do

echo ${p1[i]} ${p2[i]} $i

done

printf "\n"

mkdir -p $ID/raw_data

aws s3 --no-sign-request cp --recursive --include="*.fastq.gz" --exclude="*I1*" s3://genomeark/species/${SPECIES}/${ID}/genomic_data/10x/ ${ID}/raw_data

dw_date=`date "+%Y%m%d-%H%M%S"`
printf "Data downloaded at: $dw_date\n\n"

mkdir -p $ID/trimmed

cat ${ID}/raw_data/*R1*.fastq.gz > $ID/trimmed/fw.fastq.gz
cat ${ID}/raw_data/*R2*.fastq.gz > $ID/trimmed/rv.fastq.gz

rm ${ID}/raw_data/*

dw_date=`date "+%Y%m%d-%H%M%S"`
printf "Reads collated at: $dw_date\n\n"

/bin/proc10xG/process_10xReads.py -a -1 $ID/trimmed/fw.fastq.gz -2 $ID/trimmed/rv.fastq.gz -o $ID/trimmed/trimmed

mv $ID/trimmed/trimmed_R1_001.fastq.gz $ID/trimmed/fw.fastq.gz
mv $ID/trimmed/trimmed_R2_001.fastq.gz $ID/trimmed/rv.fastq.gz

dw_date=`date "+%Y%m%d-%H%M%S"`
printf "Barcodes trimmed at: $dw_date\n\n"

cp config.txt $ID

mkdir -p $ID/ref

cp refs/${REF}* $ID/ref

cd $ID

sed -i "s/Test/${ID}/g" config.txt

ref=$(ls ref)

sed -i "s/ref.fasta/ref\/${ref}/g" config.txt

sed -i "s/read1.fq/trimmed\/fw.fastq.gz/g" config.txt
sed -i "s/read2.fq/trimmed\/rv.fastq.gz/g" config.txt

perl ../NOVOPlasty/NOVOPlasty3.8.3.pl -c config.txt

dw_date=`date "+%Y%m%d-%H%M%S"`
printf "Completed at: $dw_date\n\n"
