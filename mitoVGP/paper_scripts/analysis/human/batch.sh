#!/bin/bash

while read bam; do

sbatch bam2fastq.sh ${bam}

done<${1}
