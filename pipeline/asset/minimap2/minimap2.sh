#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./minimap2 <ref> [num_line num_cpus]"
	echo "    Assumes we have <ref> and <input.fofn> in the same dir"
	exit -1
fi

module load minimap2/2.17
module load samtools

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

# Unless specified, use slurm array task id for input line num.
if [ -z $2 ]; then
	i=$SLURM_ARRAY_TASK_ID
else
	i=$2
fi

qry=`sed -n ${i}p input.fofn`
out=`basename $qry`
out=`echo $out | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

if [ ! -z $3 ]; then
	cpus=$3
else
	cpus=$((SLURM_CPUS_PER_TASK-1))
fi

if [ -s $out.paf ]; then
	echo "$out.paf found. Skip alignment."
else
	echo "Start aligning $qry to $ref.idx"

	echo "\
	minimap2 -x map-pb -t $cpus $ref.idx $qry > $out.paf"
	minimap2 -x map-pb -t $cpus $ref.idx $qry > $out.paf
fi

