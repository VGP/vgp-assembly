#!/bin/bash

id=$1

seq1=$2

seq2=$3

mkdir -p $id

makeblastdb -in $seq1 -parse_seqids -dbtype nucl -out $id/blast.db

blastn -outfmt "1" -query $seq2 -db $id/blast.db > $id/default_${id}.out

blastn -outfmt "6 qseqid sseqid qlen slen length pident" -query $seq2 -db $id/blast.db  | sort -k3 -n > $id/output_${id}.out

printf "qseqid\tsseqid\tqlen\tslen\tlength\tpident\tpident_noiupac\n" > $id/avg_output_${id}.out

IUPAC=$(cat $seq1 | grep -v ">" | tr -cd 'R|Y|S|W|K|M|B|D|H|V|N' | wc -c)

echo $IUPAC

awk -v IUPAC="$IUPAC" '{if (len<$4) {len+=$5; pident+=$5*$6}}END{print $1"\t"$2"\t"$3"\t"$4"\t"len"\t"(pident / len)"\t"(pident / len) + (IUPAC / len * 100)}' $id/output_${id}.out > $id/avg_output_${id}.out
