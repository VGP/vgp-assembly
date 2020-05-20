#!/bin/bash

if [[ -z $1 ]] || [[ ! -e input.fofn ]] ; then
	echo "Usage: ./_submit_bwa.sh <ref.fasta> [jid]"
	echo -e "\t<ref.fasta>: no gz!"
	echo -e "\tAssumes <input.fofn> and gaps.bed in the same path"
	echo -e "\t[jid]: set --dependency=afterok:[jid]"
	exit -1
fi

#export asset=$tools/asset

ref=$1
ref=`basename $ref`
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g'`

mkdir -p logs

if [[ ! -s $1.sa ]] && [[ -z $2 ]]; then
	cpus=4
	mem=32g
	name=index.$ref
	script=$VGP_PIPELINE/asset/bwa/index.sh
	args=$1
	walltime=1-0
	log=logs/$name.%A.log
	echo "\
	sbatch --partition=norm -D `pwd` --job-name=$name --time=$walltime --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
	sbatch --partition=norm -D `pwd` --job-name=$name --time=$walltime --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args > idx_jid
	wait_for="--dependency=afterok:"`cat idx_jid`
fi

if [ ! -z $2 ]; then
	wait_for="--dependency=afterok:$2"
fi

if [ ! -s hic.bed ]; then
	LEN=`wc -l input.fofn | awk '{print $1}'`

	if [ ! -s $LEN.bam ]; then
		cpus=32
		mem=32g
		name=hic.map.$ref
		script=$VGP_PIPELINE/asset/hic/map_hic.sh
		args="$1 input.fofn"
		log=logs/$name.%A_%a.log
	
		echo "\
		sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args"
		sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args > bwa_jid
	fi

	cpus=4
	mem=10g
	name=hic.ast.$ref
	script=$VGP_PIPELINE/asset/hic/ast.sh
	args=gaps.bed
	walltime=4:00:00
	dependency=""
	if [ -e bwa_jid ]; then
	        dependency=`cat bwa_jid`
	        dependency="--dependency=afterok:$dependency"
	fi
	log=$PWD/logs/$name.%A.log
	echo "\
	sbatch --partition=quick --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=quick --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
fi

