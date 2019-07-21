#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./minimap2 <ref>"
	echo "    Assumes we have <ref> and <input.fofn> in the same dir"
	exit -1
fi

module load minimap2/2.17
module load samtools

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

cpus=$((SLURM_CPUS_PER_TASK-1))

# Unless specified, use slurm array task id for input line num.
if [ -z $2 ]; then
	i=$SLURM_ARRAY_TASK_ID
else
	i=$2
fi

qry=`sed -n ${i}p input.fofn`

out=`basename $qry`
out=`echo $out | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

if [ -s $out.paf.gz ]; then
	echo "$out.paf.gz found. Skip alignment."
else
	echo "Start aligning $qry to $ref.idx"

	echo "\
	minimap2 -x map-pb -t $cpus $ref.idx $qry | gzip -c - > $out.read.paf.gz"
	minimap2 -x map-pb -t $cpus $ref.idx $qry | gzip -c - > $out.read.paf.gz
fi

