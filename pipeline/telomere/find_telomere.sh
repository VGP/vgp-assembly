#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./find_telomere.sh <fasta>"
	exit -1
fi

file=$1
file_name=`basename $file`

if [ ! -e $file_name ]; then
	ln -s $file
fi

file=$file_name
prefix=`echo $file | sed 's/.fasta$//g' | sed 's/.fa$//g'`

module load minimap2	# For sdust

$VGP_PIPELINE/telomere/find_telomere $file | awk '{print $1"\t"$(NF-4)"\t"$(NF-3)"\t"$(NF-2)"\t"$(NF-1)"\t"$NF}' - > $prefix.telomere
sdust $file > $prefix.sdust
java -cp $VGP_PIPELINE/telomere/telomere.jar SizeFasta $file > $prefix.lens

# Lowering threshold to 0.10 (10%) from the initial 0.40 (40%)
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereWindows $prefix.telomere 99.9 0.1 > $prefix.windows
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereBreaks $prefix.lens $prefix.sdust $prefix.telomere > $prefix.breaks


