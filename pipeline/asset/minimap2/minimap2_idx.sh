#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./minimap2_idx.sh <ref> [num_cpus]"
	echo "Assumes <ref is a fasta file"
	exit -1
fi

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g'`

if [ ! -z $2 ]; then
        cpus=$2
else
        cpus=$SLURM_CPUS_PER_TASK
fi

if [ -f $ref.idx ]; then
	echo "$ref.idx found Exit."
	exit 0
fi

echo "Start indexing $1"

module load minimap2/2.17

echo "\
minimap2 -t $cpus -x map-pb -d $ref.idx $1"
minimap2 -t $cpus -x map-pb -d $ref.idx $1
echo

