#! /bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_mashmap.sh <ref.fasta> <qry.fasta>"
	exit -1
fi

fasta1=$1
fasta2=$2

if [ ! -e `basename $fasta1` ]; then
	ln -s $fasta1
fi

if [ ! -e `basename $fasta2` ]; then
	ln -s $fasta2
fi

genome1=${fasta1/.fasta/}
genome2=${fasta2/.fasta/}

cpus=32
mem=24g
name=${genome1}_mashmap2
#script=$VGP_PIPELINE/mashmap/mashmap_small.sh
#script=$VGP_PIPELINE/mashmap/mashmap_pi90.sh
script=$VGP_PIPELINE/mashmap/mashmap.sh
ref=${genome1}
qry=${genome2}
args="$ref $qry"

mkdir -p logs
log=logs/$name.%A_%a.log

echo "" > mashmap2_jid

if [ ! -e mashmap_${genome2}_to_${genome1}/out.map ]; then 
	echo "\
	sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log -D `pwd` $script $args"
	sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log -D `pwd` $script $args > mashmap2_jid
fi

