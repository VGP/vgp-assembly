#!/bin/bash

if [ -z $1 ] ; then
	echo "Usage: ./longranger.sh <genome_id>"
	exit -1
fi

genome=$1
ref=${genome}_t1.fasta

if ! [ -e $ref ]; then
    ln -s ../../$ref
fi

ref=${ref/.fasta/}      ## if contains .fasta, remove it

if ! [ -e "refdata-$ref/genome" ] ; then
echo "No reference found"
echo "=== start indexing reference ==="
echo "\
/data/Phillippy/tools/longranger/longranger-2.2.2/longranger mkref $ref.fasta"
/data/Phillippy/tools/longranger/longranger-2.2.2/longranger mkref $ref.fasta
echo ""
fi

echo "=== start running longranger wgs ==="
/data/Phillippy/tools/longranger/longranger-2.2.2/longranger align \
--id=$genome \
--fastq=/data/rhiea/genome10k/$genome/genomic_data/10x/ \
--sample=$genome \
--reference=refdata-$ref \
--jobmode=slurm \
--localcores=32 \
--localmem=60 \
--maxjobs=500 \
--disable-ui \
--nopreflight

