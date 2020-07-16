#!/bin/bash

id=$1

seq1=$2

seq2=$3

mkdir -p $id

awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' $seq1 | sort -rnk2 > $id/contig_${id}.lengths

longest_contig=$(head -1 $id/contig_${id}.lengths | awk '{print $1}')

awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $seq1 | grep "$longest_contig" -A1 | head -2 > $id/longest_contig_${id}.fasta

makeblastdb -in $id/longest_contig_${id}.fasta -parse_seqids -dbtype nucl -out $id/blast.db

blastn -outfmt "1" -query $seq2 -db $id/blast.db > $id/default_${id}.out

blastn -outfmt "6 qseqid sseqid qlen slen length pident" -query $seq2 -db $id/blast.db  | sort -k3 -n > $id/output_${id}.out

printf "qseqid\tsseqid\tqlen\tslen\tlength\tpident\tpident_noiupac\n" > $id/avg_output_${id}.out

IUPAC=$(cat $seq1 | grep -v ">" | tr -cd 'R|Y|S|W|K|M|B|D|H|V|N' | wc -c)

echo $IUPAC

awk -v IUPAC="$IUPAC" '{if (len<$4) {len+=$5; pident+=$5*$6}}END{print $1"\t"$2"\t"$3"\t"$4"\t"len"\t"(pident / len)"\t"(pident / len) + (IUPAC / len * 100)}' $id/output_${id}.out > $id/avg_output_${id}.out
