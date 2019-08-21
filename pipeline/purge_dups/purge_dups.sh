#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: pbstats.sh <pri_asm.fasta>"
	exit -1
fi

pri_asm=$1

export PATH=$tools/purge_dups/bin:$PATH

echo "
pbcstat *.read.paf.gz "
pbcstat *.read.paf.gz  #&&

#rm *.paf.gz || echo "ERROR"

calcuts PB.stat > cutoffs 2>calcults.log
purge_dups -2 -T cutoffs -c PB.base.cov $pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log

get_seqs dups.bed $pri_asm > purged.fa 2> hap.fa 


