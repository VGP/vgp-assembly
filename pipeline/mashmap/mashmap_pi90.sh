#!/bin/bash

ref=$1.fasta
qry=$2.fasta
out=mashmap_${2}_to_${1}

mkdir -p $out

out=$out/out.map

if [ ! -e $out ]; then
	echo "\
	$tools/mashmap/mashmap -r $ref -q $qry -t $SLURM_CPUS_PER_TASK -o $out --filter_mode one-to-one --pi 90"
	$tools/mashmap/mashmap -r $ref -q $qry -t $SLURM_CPUS_PER_TASK -o $out --filter_mode one-to-one --pi 90
	echo
fi

module load gnuplot

cd mashmap_${2}_to_${1}

echo "\
$tools/mashmap/generateDotPlot png large out.map"
$tools/mashmap/generateDotPlot png large out.map
mv out.png ../mashmap_${2}_to_${1}.png

cat out.map | awk -F " " -v name=$name '{print $6"\t"$8"\t"$9+1"\t"$1"\t"$3"\t"$4+1"\t"$1":"$3"-"$4+1":"name"\t"$NF"\t"$5"\t"$7"\t"$2}' > out.map.bed
