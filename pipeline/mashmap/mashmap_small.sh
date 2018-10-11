#!/bin/bash

ref=$1.fasta
qry=$2.fasta
out=mashmap_${2}_to_${1}_2k

mkdir -p $out

out=$out/out.map

if [ ! -e $out ]; then
	echo "\
	$tools/mashmap/mashmap -r $ref -q $qry -t $SLURM_CPUS_PER_TASK -o $out -s 2000 --filter_mode one-to-one"
	$tools/mashmap/mashmap -r $ref -q $qry -t $SLURM_CPUS_PER_TASK -o $out -s 2000 --filter_mode one-to-one
	echo
fi

module load gnuplot

cd mashmap_${2}_to_${1}

echo "\
$tools/mashmap/generateDotPlot png large out.map"
$tools/mashmap/generateDotPlot png large out.map
mv out.png ../mashmap_${2}_to_${1}.png

