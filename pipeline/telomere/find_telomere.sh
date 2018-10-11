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

$VGP_PIPELINE/telomere/find_telomere $file > $prefix.telomere
sdust $file > $prefix.sdust
java $VGP_PIPELINE/telomere/SizeFasta $file > $prefix.lens

java $VGP_PIPELINE/telomere/FindTelomereWindows $prefix.telomere 99.9 > $prefix.windows
java $VGP_PIPELINE/telomere/FindTelomereBreaks $prefix.lens $prefix.sdust $prefix.telomere > $prefix.breaks

cat $prefix.windows |awk '{if ($4 > 100000 && $3-$5 > 100000) print $0}'
cat $prefix.breaks | awk '{if ($4 > 100000 && $NF-$6 > 100000) print $0}'

