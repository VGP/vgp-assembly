#!/bin/sh 

## This is the final script for removing terminal gaps and contaminant scaffolds. 

scaffs=$1
outfile=$2
contam_scaffs=$3

echo "gfastats $scaffs --exclude-bed $contam_scaffs -o $outfile --remove-terminal-gaps"
gfastats $scaffs --exclude-bed $contam_scaffs -o $outfile --remove-terminal-gaps 

