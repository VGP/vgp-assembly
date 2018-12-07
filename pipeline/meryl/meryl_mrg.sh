#!/bin/bash

k=$1
prefix=$2
filelist=$3

for file in $(grep --color=never $prefix $filelist);
do
	name=${file/.fastq.gz/}
	name=${name/.fq.gz/}
	indb=$indb" -s $name.$k"
done

mrg=$prefix.k$k

if ! [ -e $mrg.mcdat ]; then
	echo "\
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -M add $indb -o $mrg"
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -M add $indb -o $mrg

	echo "\
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -Dh -s $mrg > $mrg.hist"
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -Dh -s $mrg > $mrg.hist
fi

if ! [ -e $mrg.hist.ploidy ] ; then
	echo "\
	java -jar -Xmx512m $VGP_PIPELINE/meryl/kmerHistToPloidyDepth.jar $mrg.hist > $mrg.hist.ploidy"
	java -jar -Xmx512m $VGP_PIPELINE/meryl/kmerHistToPloidyDepth.jar $mrg.hist > $mrg.hist.ploidy
fi

x=`sed -n 2p $mrg.hist.ploidy | awk '{print $NF}'`

if ! [ -e $mrg.filt.mcdat ]; then
	echo "\
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -M greaterthan $x -s $mrg -o $mrg.filt"
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -M greaterthan $x -s $mrg -o $mrg.filt

	echo "\
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -Dc -s $mrg.filt > $mrg.filt.count"
	$tools/canu/canu-1.7/Linux-amd64/bin/meryl -Dc -s $mrg.filt > $mrg.filt.count
fi

echo "Use $mrg.hist to run genomescope."
echo
echo "k-mers are filtered with > $x to get the descriptive k-mers (See \'destict\' below)."
cat $mrg.filt.count
