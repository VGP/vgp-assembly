#!/bin/bash
set -e

echo "Usage: ./dedup.sh <alignment.bam> <threads>"

module load samtools/1.9
module load picard/2.9.2

bam=$1

mkdir -p sort_dedup
bam_prefix=`basename $bam`
sort_bam='sort_dedup/'${bam_prefix/.bam/.sort.bam}
dp_bam=${sort_bam/.bam/.dp.bam}
resort_bam=${dp_bam/.bam/.sort_n.bam}

threads=$2

#:<<'END'
echo "# Sorting by coordinates"
echo "\
samtools sort -@$threads -T $sort_bam.tmp -m2G -O bam -o $sort_bam $bam"
samtools sort -@$threads -T $sort_bam.tmp -m2G -O bam -o $sort_bam $bam
echo "\
samtools index $sort_bam"
samtools index $sort_bam
echo ""
#END

mkdir -p tmp
echo "\
java -jar -Xmx180g -Djava.io.tmpdir=$PWD/tmp /usr/local/apps/picard/2.9.2/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$sort_bam O=$dp_bam M=${dp_bam/.bam/.metrics.txt} ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024"
java -jar -Xmx180g -Djava.io.tmpdir=$PWD/tmp /usr/local/apps/picard/2.9.2/build/libs/picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=$sort_bam O=$dp_bam M=${dp_bam/.bam/.metrics.txt} ASSUME_SORT_ORDER=coordinate MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1024
echo ""
#END

echo "# Resort to name order"
echo "\
samtools sort -@$threads -n -T $resort_bam.tmp -m2000m -O bam -o $resort_bam $dp_bam"
samtools sort -@$threads -n -T $resort_bam.tmp -m2000m -O bam -o $resort_bam $dp_bam
echo ""

