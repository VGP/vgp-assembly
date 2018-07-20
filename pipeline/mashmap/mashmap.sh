#!/bin/bash

ref=$1.fasta
qry=$2.fasta
out=mashmap_${1}_to_${2}

mkdir -p $out

out=$out/$out.map

echo "\
/data/Phillippy/tools/mashmap/mashmap -r $ref -q $qry -t $SLURM_CPUS_PER_TASK -o $out"
/data/Phillippy/tools/mashmap/mashmap -r $ref -q $qry -t $SLURM_CPUS_PER_TASK -o $out
echo

module load gnuplot

echo "\
/data/Phillippy/tools/mashmap/generateDotPlot png large $out"
/data/Phillippy/tools/mashmap/generateDotPlot png large $out
