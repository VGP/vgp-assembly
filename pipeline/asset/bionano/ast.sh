#!/bin/bash

echo "Usage: ./ast.sh <input.fofn>"

input=$1
echo $input

if [[ -z $1 ]]; then
	echo "input.fofn required. Exit."
	exit -1
fi

fln=`wc -l $input | awk '{print $1}'`
echo "$fln bnx found."

if [[ ! -s bn.bed ]]; then
if [[ "$fln" -gt 2 ]]; then
	echo ">2 bionano enzymes found."
	ls bnx_*/bionano_*.bed > bed.list
	bed1=`sed -n 1p bed.list`
	for i in $(seq 2 $LEN);
	do
		bed2=`sed -n ${i}p bed.list`
		$asset/bin/union $bed1 $bed2 > bn$i.bed
		bed1=bn$i.bed
	done
	mv $bed1 bn.bed
	rm bn[1-$(($LEN-1))].bed
elif [[ "$fln" -eq 2 ]]; then
	$asset/bin/union bnx_*/bionano_*.bed > bn.bed
elif [[ "$fln" -eq 1 ]]; then
	echo "Copying bionano_*.bed to bn.bed"
	cp bnx_*/bionano_*.bed bn.bed
fi
fi

module load bedtools

platform=bn

cols=`head -n2 $platform.bed | tail -n1 | awk '{print NF}'`
if [[ "$cols" -gt 3 ]]; then
	awk '{print $1"\t"$(NF-1)"\t"$NF}' $platform.bed > $platform.format.bed
else
	cp $platform.bed $platform.format.bed
fi

bedtools subtract -a ../asm.bed -b $platform.format.bed | bedtools merge -d 100 -i - > $platform.low_high.bed
bedtools subtract -a $platform.low_high.bed -b ../asm.ends.bed -A > $platform.low_high.trim1k.bed
bedtools subtract -a ../asm.bed -b $platform.low_high.bed > $platform.support.bed
