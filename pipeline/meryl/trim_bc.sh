#!/bin/bash

if [ -z $SLURM_ARRAY_TASK_ID ]; then
	i=$2
else
	i=$SLURM_ARRAY_TASK_ID
fi

out_prefix=$1

r1=`sed -n ${i}p input.fofn`
r2=${r1/R1_001/R2_001}

if ! [[ -e $out_prefix.r1.$i.fastq && -e $out_prefix.$i.name ]] ; then
	if [[ "$r1" =~ \.gz$ ]]; then
		gunzip $r1
		r1=${r1/.gz/}
	fi

	echo "scaff_BC-reads-1 $r1 $out_prefix.r1.$i.fastq read.$i.name"
	$tools/scaff10x/Scaff10X/src/scaff-bin/scaff_BC-reads-1 $r1 $out_prefix.r1.$i.fastq $out_prefix.$i.name
	echo ""
fi

if ! [ -e $out_prefix.r2.$i.fastq ]; then
        if ! [[ "$r2" =~ \.gz$ ]] && [[ -e $r2.gz ]]; then
                gunzip $r2
        elif [[ "$r2" =~ \.gz$ ]]; then
                gunzip $r2
		r2=${r2/.gz/}
	fi
	echo "scaff_BC-reads-2 $out_prefix.$i.name $r2 $out_prefix.r2.$i.fastq"
	$tools/scaff10x/Scaff10X/src/scaff-bin/scaff_BC-reads-2 $out_prefix.$i.name $r2 $out_prefix.r2.$i.fastq
fi
#rm $r1
#rm $r2
#rm $out_prefix.$i.name
