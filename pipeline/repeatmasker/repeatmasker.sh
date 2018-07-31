#!/bin/bash

module load repeatmasker
### repeatmasker/4.0.7

species=vertebrates

par=$(($SLURM_CPUS_PER_TASK / 2))
if [ $par -lt 2 ]
then
    par=2
fi


fa=asm.fasta

if ! [ -e $fa ]; then
	echo "No $fa. Exit."
	exit -1
fi

echo "\
RepeatMasker -pa $par -q -species=$species -xm -dir=out $fa"
RepeatMasker -pa $par -q -species=$species -xm -dir=out $fa
