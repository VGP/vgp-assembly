#!/bin/bash

species=$1

id=$2

mkdir -p $id

aws s3 cp --exclude "*" --include "mtDNA_*.fasta" --recursive s3://genomeark/species/${species}/${id}/assembly_MT_rockefeller/intermediates/reference/ $id

GSIZE=$(awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $id/mtDNA_*.fasta | awk -F'\t' '{print $2}')

seqtk seq -a ${id}.fastq* > $id/${id}.fasta

makeblastdb -in $id/${id}.fasta -parse_seqids -dbtype nucl -out $id/blast.db

blastn -dust no -outfmt "6 sseqid slen qcovs" -query $1 -db blast.db  | sort -k2 -nr | uniq |awk -v gsize="${GSIZE}" '{printf $0 "\t" gsize/$2*$3 "\n"}' | awk -v len="${4}"  -v per="${3}" '{if($4>per && $2<len) printf $0 "\n"}' > reads${3}_${4}.out
