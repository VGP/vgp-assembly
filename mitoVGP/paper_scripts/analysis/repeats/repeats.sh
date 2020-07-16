#!/bin/bash

mkdir -p $1

cd $1

cat ../../../../VGP_dataset/$1* ../NOVOplasty/$1/Circular* > kmers.fasta

/bin/windowmasker -in kmers.fasta -mk_counts -out kmer.counts

/bin/windowmasker -in ../../../../VGP_dataset/$1* -ustat kmer.counts -out $1.report

/bin/windowmasker -in ../NOVOplasty/$1/Circular* -ustat kmer.counts -out NOVOplasty.report

awk -F' ' '{diff=$3 - $1; print diff; diffs=diffs+diff}END{print diffs}' $1.report > $1.diffs

awk -F' ' '{diff=$3 - $1; print diff; diffs=diffs+diff}END{print diffs}' NOVOplasty.report > NOVOplasty.diffs

VGP=$(tail -1 $1.diffs)
NOV=$(tail -1 NOVOplasty.diffs)

printf "%s\t%s\t%s\t%s\n" $1 $VGP $NOV $(( $VGP - $NOV )) > summary.txt
