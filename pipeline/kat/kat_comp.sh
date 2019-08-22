#!/bin/bash

source /data/Phillippy/tools/conda/etc/profile.d/conda.sh
conda activate kat

assembly=$1

threads=$SLURM_CPUS_PER_TASK

out=`echo $assembly | awk -F "/" '{print $NF}'`
out=${out/.fasta/}
mkdir -p $out
out=$out/kat_comp

options="-H 10000000000 -I 10000000000 -m 21 -h"

echo "\
kat comp -o $out -t $threads $options 'R1.fastq.gz R2.fastq.gz' $assembly"
kat comp -o $out -t $threads $options 'R1.fastq.gz R2.fastq.gz' $assembly

