#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ast.sh <out_prefix> [mean_cov]"
	exit -1
fi

prefix=$1
mean_cov=$2
postfix=""

if [[ ! -z $mean_cov ]]; then
    mean_cov=`echo $mean_cov | awk '{printf "%.0f\n", $1*2.5}'`
    mean_cov="-M $mean_cov"
    postfix="_M"
fi

#export asset=$tools/asset
pafs=`ls *.paf`

bed=pb${postfix}.bed

if [[ ! -s $bed ]]; then
echo "
	$asset/bin/ast_pb $mean_cov $pafs > $bed"
	$asset/bin/ast_pb ${mean_cov} $pafs > $bed && echo "Asset Finished!" || exit -1
fi

echo "clean up"
#rm *.paf

module load bedtools

bedtools subtract -a ../asm.bed -b $bed | bedtools merge -d 100 -i - > pb.low_high.bed

# Remove ends: remove low coverage region around ~1kb
bedtools subtract -a pb.low_high.bed -b ../asm.ends.bed -A > pb.low_high.trim1k.bed

# Get support: for getting reliable blocks. Will trim at the end when getting acc.reliable.bed
bedtools subtract -a ../asm.bed -b pb.low_high.bed > pb.support.bed

