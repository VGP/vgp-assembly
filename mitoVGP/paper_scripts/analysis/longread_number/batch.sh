#!/bin/bash

while IFS=$'\t' read -r species id
do

mkdir -p $id

aws s3 cp --exclude "*" --include "*_MT_extracted_reads/${id}.fastq*" --recursive s3://genomeark/species/${species}/${id}/assembly_MT_rockefeller/intermediates/ $id

aws s3 cp --exclude "*" --include "mtDNA_*.fasta" --recursive s3://genomeark/species/${species}/${id}/assembly_MT_rockefeller/intermediates/reference/ $id

mv $id/*_MT_extracted_reads/* $id
rm -r $id/*_MT_extracted_reads/

GSIZE=$(awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $id/mtDNA_*.fasta | awk -F'\t' '{print $2}')

seqtk seq -a $id/${id}.fastq* > $id/${id}.fasta

makeblastdb -in $id/${id}.fasta -parse_seqids -dbtype nucl -out $id/blast.db

blastn -outfmt "6 sseqid slen qcovs" -query $id/mtDNA_*.fasta -db $id/blast.db  | sort -k2 -nr | uniq |awk -v gsize="${GSIZE}" '{printf $0 "\t" gsize/$2*$3 "\n"}' | awk  -v per="${1}" '{if($4>per) printf $0 "\n"}' > $id/reads_${id}.out

done<species.ls
