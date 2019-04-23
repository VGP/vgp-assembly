#!/bin/bash

asm=$1

if [ -z $asm ]; then
	echo "Usage: ./trimNs.sh <asm.fasta>"
	exit -1
fi

echo "Start running trim-N process on $asm"

name=`echo $asm | sed 's/.fasta$//g' | sed 's/.fa$//g'`

module load python/3.5

echo "
python3 $VGP_PIPELINE/bionano/trimNs/remove_fake_cut_sites_DNAnexus.py $asm $name.rmFake.fasta remove_fake_cuts.log
"
python3 $VGP_PIPELINE/bionano/trimNs/remove_fake_cut_sites_DNAnexus.py $asm $name.rmFake.fasta remove_fake_cuts.log

echo "
python3 $VGP_PIPELINE/bionano/trimNs/trim_Ns_DNAnexus.py $name.rmFake.fasta $name.rmFake.list
"
python3 $VGP_PIPELINE/bionano/trimNs/trim_Ns_DNAnexus.py $name.rmFake.fasta $name.rmFake.list

echo "
python3 $VGP_PIPELINE/bionano/trimNs/clip_regions_DNAnexus.py $name.rmFake.fasta $name.rmFake.list $name.trimmed.fasta
"
python3 $VGP_PIPELINE/bionano/trimNs/clip_regions_DNAnexus.py $name.rmFake.fasta $name.rmFake.list $name.trimmed.fasta

