#! /bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_mashmap.sh <job-name> <ref> <qry>"
	exit -1
fi

cpus=32
mem=12g
name=${1}_mashmap
script=/data/Phillippy/tools/vgp-assembly/git/vgp-assembly/pipeline/mashmap/mashmap.sh
ref=$2
qry=$3

ref=${ref/.fasta/}
qry=${qry/.fasta/}
args="$ref $qry"

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

echo "\
sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args
