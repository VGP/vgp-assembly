#!/bin/bash

#this script (mtDNApipe.sh) is used to retrieve mitochondrial-like sequences from the raw Pacbio data 
#generated in the framework of the Vertebrate Genomes Project and assemble them using Canu.

#it requires the following software (and their dependencies) installed:
#aws-cli/1.16.101, blasr/5.3.2-06c9543, bedtools2/2.18 (bamToFastq), Canu/1.8, blastn/2.7.1+

#sequence retrieval is based on a search by similarity using BLASR alignment.
#all Pacbio raw data files are individually aligned to a reference genome provided by the user.
#the reference genome can be from the same species if available, or from a
#closely-to-distantly related species.

#the approach is similar to that of Organelle_PBA described in:
#Soorni et al. BMC Genomics (2017) DOI 10.1186/s12864-016-3412-9

#in the second steps reads are used to generate assemblies using Canu assembler
#usually with default parameters (for rGopEvg1 and mRhiFer1 minOverlapLength=300 
#correctedErrorRate=0.105 were used, respectively).

#the reference genome provided by the user is then blasted to the contigs generated 
#by Canu to identify the putative mitocontig.

#required inputs are:
#1) the species name (e.g. Calypte_anna)
#2) the VGP species ID (e.g. bCalAnn1)
#3) the reference sequence fasta file
#4) the putative mitogenome size (potentially, that of the reference genome). It does not
#need to be precise. Accepts Canu formatting.
#5) the number of threads

set -e

#set variables species abbreviation size
SPECIES=$1
ABBR=$2
REF=$3
SIZE=$4
NPROC=$5

#define working directory
W_URL=${SPECIES}/assembly_MT/intermediates

if ! [[ -e "${W_URL}" ]]; then

mkdir -p ${W_URL}

fi

#copy the user-provided reference mitogenome to the reference folder
if ! [[ -e "${W_URL}/reference" ]]; then

mkdir ${W_URL}/reference
cp ${REF} ${W_URL}/reference/${REF%.*}.fasta

fi

if ! [[ -e "${W_URL}/log" ]]; then

mkdir ${W_URL}/log

fi

#record Pacbio raw data files available in the cloud at the time of the analysis
dw_date=`date "+%Y%m%d-%H%M%S"`; aws s3 --no-sign-request ls s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/pacbio/ | grep -oP "m.*.subreads.bam" | uniq > ${W_URL}/log/file_list_$dw_date.txt

if ! [[ -e "${W_URL}/pacbio_bam" ]]; then

mkdir ${W_URL}/pacbio_bam

fi

#for each Pacbio raw data file do
while read p; do

if ! [[ $p == *scraps* ]] && ! [[ $p == *.pbi ]] && [[ $p == *.bam ]] && ! [[ -e "${W_URL}/pacbio_bam/aligned_${p}" ]]; then

#download
aws s3 --no-sign-request cp s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/pacbio/$p ${W_URL}/
#align
blasr --bestn 1 ${W_URL}/$p ${W_URL}/reference/${REF%.*}.fasta --bam --out ${W_URL}/pacbio_bam/aligned_$p --nproc ${NPROC}
#remove
rm ${W_URL}/$p

fi

done <${W_URL}/log/file_list_$dw_date.txt

#organize the files

#convert to fastq
if ! [[ -e "${W_URL}/pacbio_MT_extracted_reads" ]]; then

mkdir ${W_URL}/pacbio_MT_extracted_reads

for f in ${W_URL}/pacbio_bam/aligned_*.bam; do filename=$(basename -- "$f"); filename="${filename%.*}"; bamToFastq -i $f -fq "${W_URL}/pacbio_MT_extracted_reads/${filename}.fq"; done

#merge into a single read file
if ! [[ -e "${W_URL}/pacbio_MT_extracted_reads/${ABBR}.fastq" ]]; then

cat ${W_URL}/pacbio_MT_extracted_reads/*.fq > ${W_URL}/pacbio_MT_extracted_reads/${ABBR}.fastq

fi

rm ${W_URL}/pacbio_MT_extracted_reads/*.fq

fi

#assemble mtDNA reads with canu
if ! [[ -e "${W_URL}/canu" ]]; then

canu \
 -p ${ABBR} -d ${W_URL}/canu \
 genomeSize=${SIZE} \
 -pacbio-raw ${W_URL}/pacbio_MT_extracted_reads/${ABBR}.fastq

fi

if ! [[ -e "${W_URL}/blast" ]]; then

mkdir -p ${W_URL}/blast

fi

#build blast db
makeblastdb -in ${W_URL}/canu/${ABBR}.contigs.fasta -parse_seqids -dbtype nucl -out ${W_URL}/blast/${ABBR}.db

#search the putative mitocontig using blastn
blastn -query ${W_URL}/reference/${REF%.*}.fasta -db ${W_URL}/blast/${ABBR}.db -out ${W_URL}/blast/${ABBR}.out
blastn -query ${W_URL}/reference/${REF%.*}.fasta -db ${W_URL}/blast/${ABBR}.db -outfmt 6 -out ${W_URL}/blast/${ABBR}.tb
sed -i "1iquery_acc.ver\tsubject_acc.ver\t%_identity\talignment_length\tmismatches\tgap_opens\tq.start\tq.end\ts.start\ts.end\tevalue\tbitscore" ${W_URL}/blast/${ABBR}.tb
cat ${W_URL}/blast/${ABBR}.tb | column -t > ${W_URL}/blast/${ABBR}_results.txt