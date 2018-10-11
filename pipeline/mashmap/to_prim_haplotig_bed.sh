#!/bin/bash

genome=$1

if [ -z $genome ]; then
	echo "Usage: ./to_purge_haplotigs_bed.sh <genome_id>"
	exit -1
fi

map1=mashmap_${genome}_curated_to_${genome}_v1.p/out.map
map2=mashmap_${genome}_curated.haplotigs_to_${genome}_v1.p/out.map

# Transform to bedpe

name="PRIM"
cat $map1 | awk -F " " -v name=$name '{print $6"\t"$8"\t"$9+1"\t"$1"\t"$3"\t"$4+1"\t"$1":"$3"-"$4+1":"name"\t"$NF"\t"$5"\t"$7"\t"$2}' > ${genome}_v1.p_mashmap_purge_haplotigs.bed
echo

name="HAPLOTIG"
cat $map2 | awk -F " " -v name=$name '{print $6"\t"$8"\t"$9+1"\t"$1"\t"$3"\t"$4+1"\t"$1":"$3"-"$4+1":"name"\t"$NF"\t"$5"\t"$7"\t"$2}' >> ${genome}_v1.p_mashmap_purge_haplotigs.bed
echo

echo "Done!"
