#!/bin/bash

set -e

SPECIES=$1
ABBR=$2
CONTIG=$3

cd ${SPECIES}

mkdir -p ${ABBR}_arrow

awk '$2 == "'${CONTIG}'" {print $1}' ${ABBR}_canu/${ABBR}.contigs.layout.readToTig > ${ABBR}_arrow/${ABBR}_contig${CONTIG}_IDs.txt

gunzip -c ${ABBR}_canu/${ABBR}.trimmedReads.fasta.gz > ${ABBR}_arrow/${ABBR}.trimmedReads.fasta

READS=$(grep -c ">" ${ABBR}_arrow/${ABBR}.trimmedReads.fasta)

echo "$READS reads were trimmed by Canu"

while read IDs; do grep -A 1 "id=\b$IDs\b" ${ABBR}_arrow/${ABBR}.trimmedReads.fasta; done < ${ABBR}_arrow/${ABBR}_contig${CONTIG}_IDs.txt > ${ABBR}_arrow/${ABBR}.trimmedReads.contig${CONTIG}.fasta

grep -o -P '(?<=>)(\S*)' ${ABBR}_arrow/${ABBR}.trimmedReads.contig${CONTIG}.fasta | uniq -u > ${ABBR}_arrow/${ABBR}_contig${CONTIG}_names.txt

READS_AS_U=$(grep -c "m" ${ABBR}_arrow/${ABBR}_contig${CONTIG}_names.txt)

echo "of which $READS_AS_U were used in the assembly of contig ${CONTIG}"

sed ':a;N;/>/!s/\n//;ta;P;D' ${ABBR}_canu/${ABBR}.contigs.fasta | grep -A1 ">tig0*${CONTIG} " > ${ABBR}_arrow/${ABBR}.contig${CONTIG}.fasta

for i in aligned_pb_raw_reads/*.subreads.bam; do pbindex $i; done

dataset create --type SubreadSet --name ${ABBR} ${ABBR}_arrow/${ABBR}_read_set.xml aligned_pb_raw_reads/*.subreads.bam
pbalign --maxHits 1 ${ABBR}_arrow/${ABBR}_read_set.xml ${ABBR}_arrow/${ABBR}.contig${CONTIG}.fasta ${ABBR}_arrow/${ABBR}.realigned_raw_reads.bam --nproc 16

cp aligned_pb_raw_reads/*.subreads.bam.pbi ${ABBR}_arrow/

samtools view -H ${ABBR}_arrow/${ABBR}.realigned_raw_reads.bam | sed "s/@RG/@RG\tSM:unknown/g" | sed "s/SO:UNKNOWN/SO:unknown/g" | samtools reheader - ${ABBR}_arrow/${ABBR}.realigned_raw_reads.bam > ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh.bam
java -jar /home/gformenti/bin/picard.jar FilterSamReads I=${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh.bam O=${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}.bam READ_LIST_FILE=${ABBR}_arrow/${ABBR}_contig${CONTIG}_names.txt FILTER=includeReadList VALIDATION_STRINGENCY=STRICT

samtools faidx ${ABBR}_arrow/${ABBR}.contig${CONTIG}.fasta

pbindex ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}.bam

variantCaller -j8 ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}.bam -r ${ABBR}_arrow/${ABBR}.contig${CONTIG}.fasta -o ${ABBR}_arrow/${ABBR}.contig${CONTIG}_arrow.fasta --algorithm=arrow --numWorkers 16

pbalign --maxHits 1 ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}.bam ${ABBR}_arrow/${ABBR}.contig${CONTIG}_arrow.fasta ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}_pl.bam --nproc 16

samtools faidx ${ABBR}_arrow/${ABBR}.contig${CONTIG}_arrow.fasta

pbindex  ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}_pl.bam

variantCaller -j8 ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}_pl.bam -r ${ABBR}_arrow/${ABBR}.contig${CONTIG}_arrow.fasta -o ${ABBR}_arrow/${ABBR}.contig${CONTIG}_arrow2.fasta --algorithm=arrow --numWorkers 16

pbalign --maxHits 1 ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}_pl.bam ${ABBR}_arrow/${ABBR}.contig${CONTIG}_arrow2.fasta ${ABBR}_arrow/${ABBR}.realigned_raw_reads_rh_contig${CONTIG}_pl2.bam --nproc 16

mkdir ${ABBR}_arrow/subreads_pbi/
mv ${ABBR}_arrow/*subreads.bam.pbi ${ABBR}_arrow/subreads_pbi/