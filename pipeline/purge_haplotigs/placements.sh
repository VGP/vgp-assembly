#!/bin/bash

module load purge_haplotigs
module load mummer

genome=$1

threads=$SLURM_CPUS_PER_TASK

echo "\
purge_haplotigs  ncbiplace -p ${genome}_v1.p.fasta -h ${genome}_v1.h.fasta -o ${genome}_v1.placements.tsv -t $threads -c 40"
purge_haplotigs  ncbiplace -p ${genome}_v1.p.fasta -h ${genome}_v1.h.fasta -o ${genome}_v1.placements.tsv -t $threads -c 40
