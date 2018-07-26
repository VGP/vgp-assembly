#!/bin/bash

if [ -z $1 ]; then
    echo "Usage: ./_submit_qv.sh <sample_id> [jobid_to_set_dependency]"
    exit -1
fi


sample=$1

if ! [ -e aligned.bam ]; then
        ln -s $sample/outs/possorted_bam.bam aligned.bam
        ln -s $sample/outs/possorted_bam.bam.bai aligned.bam.bai
fi

mkdir -p logs


cpus=4
mem=4g
name=$sample.genomecov
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/qv/genomecov.sh
args=$sample
walltime=1-0
log=logs/$name.%A_%a.log

if ! [ -z $2 ]; then
        wait_for="--dependency=afterok:$2"
fi
echo "\
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > genomecov_jid

