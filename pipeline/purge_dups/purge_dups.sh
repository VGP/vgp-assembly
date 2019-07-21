#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: pbstats.sh <pri_asm.fasta>"
	exit -1
fi

pri_asm=$1

export PATH=$tools/purge_dups/bin:$PATH

echo "
pbcstat *.paf.gz "
pbcstat *.paf.gz  #&&

#rm *.paf.gz || echo "ERROR"

echo "\
split_fa $pri_asm > $pri_asm.split"
split_fa $pri_asm > $pri_asm.split

minimap2 -xasm5 -DP $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz

purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

get_seqs dups.bed $pri_asm > purged.fa 2> hap.fa 


