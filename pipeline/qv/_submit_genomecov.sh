#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./_submit_genomecov.sh <pos_sorted.bam> [jobid_to_set_dependency]"
    exit -1
fi


bam=$1
sample=${bam/.bam/}

mkdir -p logs


cpus=4
mem=4g
name=$sample.genomecov
script=$VGP_PIPELINE/qv/genomecov.sh
args="$sample $sample.bam"
walltime=2-0
log=logs/$name.%A_%a.log

if ! [ -z $2 ]; then
        wait_for="--dependency=afterok:$2"
fi
echo "\
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > genomecov_jid

