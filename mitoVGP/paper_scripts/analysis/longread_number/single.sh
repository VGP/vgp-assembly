#!/bin/bash

GSIZE=$(awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $1 | awk -F'\t' '{print $2}')

seqtk seq -a $2 > ${2%.*}.fasta

makeblastdb -in ${2%.*}.fasta -parse_seqids -dbtype nucl -out blast.db

blastn -dust no -outfmt "6 sseqid slen qcovs" -query $1 -db blast.db  | sort -k2 -nr | uniq |awk -v gsize="${GSIZE}" '{printf $0 "\t" gsize/$2*$3 "\n"}' | awk -v len="${4}"  -v per="${3}" '{if($4>per && $2<len) printf $0 "\n"}' > reads${3}_${4}.out
