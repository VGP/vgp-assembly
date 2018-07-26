map=$1
name=$2

# Transform to bedpe
cat $map | awk -F " " -v name=$name '{print $6"\t"$8"\t"$9+1"\t"$1"\t"$3"\t"$4+1"\t"$1":"$3"-"$4+1":"name"\t"$NF"\t"$5"\t"$7"\t"$2}' > $map.bed

