#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./minimap2_self.sh <ref.fasta>"
	exit -1
fi

pri_asm=$1
cpus=$((SLURM_CPUS_PER_TASK-1))

PATH=$tools/purge_dups/bin:$PATH
module load minimap2

echo "\
split_fa $pri_asm > $pri_asm.split"
split_fa $pri_asm > $pri_asm.split

echo "\
minimap2 -xasm5 -DP -t $cpus $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz"
minimap2 -xasm5 -DP -t $cpus $pri_asm.split $pri_asm.split | gzip -c - > $pri_asm.split.self.paf.gz


