#!/bin/bash

if [[ -z $1 ]]; then
        echo "Usage: reliable_block.sh <genome> <genomesize.map> [chromosome_assignments.csv]"
        echo
        echo "Run after asset.sh has fully finished its submitted jobs."
        echo "Required files: asm.bed, asm.ends.bed, gaps.bed, and */*.support.bed"
        echo
        echo "    <genome>: top level folder of asset"
        echo "    <genomesize.map>: genome <tab> genome_size"
        echo "    [chromosome_assignments.csv]: Scaffold <tab> Assigned_Chromosome <tab> Localized(y/n)"
        exit 0
fi

genome=$1
map=$2 # genome <tab> genome_size
csv=$3

map=`realpath $map`

ln -s $map
map=`basename $map`


#:<<'END'
echo "Run asset"
$asset/bin/acc gaps.bed */*.support.bed > acc.bed 2> acc.log
echo

echo "Merge to get reliable blocks"
awk '$4>1' acc.bed | bedtools merge -i - > acc.gt2.mrg.bed

module load bedtools

bedtools subtract -a asm.bed -b acc.gt2.mrg.bed | bedtools merge -d 100 -i - > low_support.bed
bedtools subtract -a asm.bed -b low_support.bed > reliable.bed
bedtools subtract -A -a low_support.bed -b asm.ends.bed > low_support.trim1k.bed
#END

genome_size=`grep $genome $map | awk '{print $2}'`

echo "Get reliable block NG50 for $genome using genome size $genome_size"
awk '{print ($3-$2)}' reliable.bed | java -jar -Xmx1g $VGP_PIPELINE/stats/lenCalcNGStats.jar - $genome_size 1 | \
  awk -v name=$genome '{print name"\t"$2"\t"$3"\t"$4"\t"$14}' | head -n2 | tail -n1 >> ../reliable.low_high_100mrg.stats

if [[ -z $csv ]]; then
    echo "No .csv provided. Look if exists in assembly_curated/"
    csv=`ls /data/rhiea/genome10k/assembly_curated/$genome.pri.cur.*.chromosomes.csv`
    if [[ -z $csv ]] ; then
        echo "Attempt failed. Exit."
        exit 0;
    fi
fi

# $genome.chr.list : Assigned scaffolds
awk -F "," '{print $1}' $csv  > $genome.assigned

echo " # Reliable block NG stats for chr assigned scaffolds"
java -jar -Xmx1g $VGP_PIPELINE/utils/txtContains.jar reliable.bed $genome.assigned 1 > reliable.chr.bed
awk '{print ($3-$2)}' reliable.chr.bed | java -jar -Xmx1g $VGP_PIPELINE/stats/lenCalcNGStats.jar - $genome_size 1 | \
  awk -v name=$genome '{print name"\t"$2"\t"$3"\t"$4"\t"$14}' | head -n2 | tail -n1 >> ../reliable.low_high_100mrg.assigned.stats


