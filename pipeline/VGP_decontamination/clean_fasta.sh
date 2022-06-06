#!/bin/sh 

## This is the final script for removing terminal gaps and contaminant scaffolds. 

scaffs=$1
outfile=$2
contam_scaffs=$3
mito_scaffs=$4
newdir=$5

cat $mito_scaffs $contam_scaffs | tr [A-Z] [a-z] > $newdir/full_contam_list.txt

echo "gfastats $scaffs --exclude-bed $contam_scaffs -o $outfile --remove-terminal-gaps"
gfastats $scaffs --exclude-bed $newdir/full_contam_list.txt -o $outfile --remove-terminal-gaps 
