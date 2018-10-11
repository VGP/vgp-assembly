#! /bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_mashmap.sh <job-name>"
	exit -1
fi

genome=$1

if [ -e ../purge_haplotigs/${genome}_curated.fasta ]; then
	ln -s ../purge_haplotigs/${genome}_curated.fasta 
	ln -s ../purge_haplotigs/${genome}_curated.haplotigs.fasta
fi
if [ -e ../../${genome}_v1.p.fasta ]; then
	ln -s ../../${genome}_v1.p.fasta
else
	echo "No ${genome}_v1.p.fasta found. Exit -1."
	exit -1
fi

cpus=32
mem=24g
name=${genome}_mashmap2
script=$VGP_PIPELINE/mashmap/mashmap.sh
ref=${genome}_v1.p
qry=${genome}_curated
args="$ref $qry"

mkdir -p logs
log=logs/$name.%A_%a.log

echo "" > mashmap2_jid

if [ ! -e mashmap_${genome}_curated_to_${genome}_v1.p/out.map ]; then 
	echo "\
	sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log -D `pwd` $script $args"
	sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log -D `pwd` $script $args > mashmap2_jid
fi

qry=${genome}_curated.haplotigs
args="$ref $qry"
log=logs/$name.%A_%a.log

if [ ! -e mashmap_${genome}_curated.haplotigs_to_${genome}_v1.p/out.map ]; then
	echo "\
	sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log -D `pwd` $script $args"
	sbatch --partition=quick --job-name=$name --time=4:00:00 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log -D `pwd` $script $args >> mashmap2_jid
fi
jid=`cat mashmap2_jid | tr '\n' ',' | sed 's/,$//g'`
wait_for="--dependency=afterok:$jid"

cpus=4
mem=8g
name=${genome}_to_bed
script=$VGP_PIPELINE/mashmap/to_prim_haplotig_bed.sh
args=$genome
log=logs/$name.%A_%a.log

echo "\
sbatch --partition=quick --job-name=$name --time=30:00 --cpus-per-task=$cpus --mem=$mem $wait_for --error=$log --output=$log -D `pwd` $script $args"
sbatch --partition=quick --job-name=$name --time=30:00 --cpus-per-task=$cpus --mem=$mem $wait_for --error=$log --output=$log -D `pwd` $script $args 
