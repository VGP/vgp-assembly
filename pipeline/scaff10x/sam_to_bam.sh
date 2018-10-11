#!/bin/bash

module load samtools

cpus=$SLURM_CPUS_PER_TASK
samtools view -hb -@$cpus align.sam > align.bam

