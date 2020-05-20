#!/bin/bash

if [[ "$#" -lt 2 ]]; then
    echo "Usage: ./nucmer.sh <ref.fa> <qry.fa>"
    exit 0
fi

ref=$1  # ref.fasta
qry=$2  # qry.fasta

ref_prefix=`basename $ref`
ref_prefix=`echo $ref_prefix | sed 's/.fasta$//g' | sed 's/.fa$//g'`
qry_prefix=`basename $qry`
qry_prefix=`echo $qry_prefix | sed 's/.fasta$//g' | sed 's/.fa$//g'`

prefix=nucmer.${qry_prefix}_to_${ref_prefix}.l100.c500

mkdir -p $prefix
out=$prefix/$prefix

module load mummer/3.23

echo "
nucmer -p $out --maxmatch -l 100 -c 500 $ref $qry"
nucmer -p $out --maxmatch -l 100 -c 500 $ref $qry

echo "
dnadiff -p $out.dnadiff -d $out.delta"
dnadiff -p $out.dnadiff -d $out.delta

echo "
mummerplot -p $out.dnadiff.plot.1 --layout --fat -t png $out.dnadiff.1delta"
mummerplot -p $out.dnadiff.plot.1 --layout --fat -t png $out.dnadiff.1delta


