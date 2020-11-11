#! /bin/bash

echo "Usage: ./_submit_bwa.sh <ref.fasta> <fastq.map> <out>"

ref=$1
fastq_map=$2
out=$3

if [[ "$#" -lt 3 ]]; then
	echo "Not enough args given. Exit."
	exit -1
fi

mkdir -p logs

if [ ! -s $ref.fai ]; then
	cpus=4
	mem=10g
	name=$out.index
	script=$VGP_PIPELINE/bwa/bwa_index.sh
	args="$ref"
	partition=quick
	walltime=4:00:00
	path=`pwd`
	log=logs/$name.%A.log

	echo "\
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
	sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > index.jid
	index=`cat index.jid`
	echo $index
	extra="--dependency=$index"
fi

cpus=24
mem=30g
name=$out.bwa
script=$VGP_PIPELINE/bwa/bwa.sh
args="$ref $fastq_map $out"
partition=norm
walltime=1-0
path=`pwd`
log=logs/$name.%A_%a.log

LEN=`wc -l $fastq_map | awk '{print $1}'`
extra="--array=1-$LEN $extra"

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args | awk '{print $NF}' > bwa.jid
wait_for="--dependency=afterok:"`cat bwa.jid`

if [[ -e bam.list ]]; then
	rm bam.list
fi
for i in $(seq 1 $LEN)
do
	echo "$out.$i.bam" >> bam.list
done

cpus=24
mem=48g
name=$out.merge
script=$VGP_PIPELINE/bwa/merge.sh
args="$out"
partition=norm
walltime=1-0
path=`pwd`
log=logs/$name.%A.log
extra=$wait_for

echo "\
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args"
sbatch -J $name --mem=$mem --partition=$partition --cpus-per-task=$cpus -D $path $extra --time=$walltime --error=$log --output=$log $script $args > mrg.jid


sample=$out

if ! [ -e aligned.bam ]; then
	ln -s $sample.bam aligned.bam
	ln -s $sample.bam.bai aligned.bam.bai
fi

mkdir -p refdata-asm/fasta
cd refdata-asm/fasta
ln -s ../../$ref genome.fa
ln -s ../../$ref.fai genome.fa.fai
cd ../../

cpus=4
mem=12g
name=$sample.freebayes
script=$VGP_PIPELINE/freebayes-polish/freebayes_v1.3.sh
args=$sample
walltime=2-0
log=logs/$name.%A_%a.log

mkdir -p bcf

wait_for="--dependency=afterok:`cat mrg.jid`"

echo "\
sbatch --partition=norm --array=1-100 $wait_for -D $PWD --cpus-per-task=$cpus -J $name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --array=1-100 $wait_for -D $PWD --cpus-per-task=$cpus -J $name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > fb.jid
wait_for="--dependency=afterok:`cat fb.jid`"

cpus=2
mem=4g
name=$sample.consensus
script=$VGP_PIPELINE/freebayes-polish/consensus.sh
args=$sample
walltime=2-0
log=logs/$name.%A_%a.log

echo "\
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > cns.jid
wait_for="--dependency=afterok:`cat cns.jid`"

cpus=4
mem=4g
name=$sample.genomecov
script=$VGP_PIPELINE/qv/genomecov.sh
args=$sample
walltime=3-0
log=logs/$name.%A_%a.log

echo "\
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > qv.jid
