#!/bin/bash

if [[ -z $1 ]] || [[ ! -e input.fofn ]] ; then
	echo "Usage: ./_submit_minimap2.sh <ref.fasta> [genome_cov]"
	echo -e "\tAssumes <input.fofn> in the same path"
	echo -e "\t[genome_cov]: set ast_pb -M to genome_cov*2.5"
	exit -1
fi

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' `
mean_cov=$2

mkdir -p logs

if [ ! -e $ref.idx ]; then
	cpus=16
	mem=32g
	name=minimap2.index.$ref
	script=$VGP_PIPELINE/asset/minimap2/minimap2_idx.sh
	args=$1
	log=logs/$name.%A.log
	echo "\
	sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
	sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args > idx_jid
	wait_for="--dependency=afterok:"`cat idx_jid`
fi

#if [ ! -s pb.bed ]; then
	LEN=`wc -l input.fofn | awk '{print $1}'`
	PAFs=`ls *.paf 2> /dev/null | wc -l | awk '{print $1}'`
	echo "$PAFs found."
	if [[ $LEN -gt $PAFs ]]; then
		cpus=16
		mem=32g
		name=minimap2.map.$ref
		script=$VGP_PIPELINE/asset/minimap2/minimap2.sh
		args=$1
		log=logs/$name.%A_%a.log
	
		echo "\
		sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args"
		sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args > minimap2_jid
	fi
#fi
	cpus=4
	mem=12g
	name=minimap2.ast.$ref
	script=$VGP_PIPELINE/asset/minimap2/ast.sh
	args="$ref $mean_cov"
	walltime=1:00:00
	dependency=""
	if [ -e minimap2_jid ]; then
	        dependency=`cat minimap2_jid`
	        dependency="--dependency=afterok:$dependency"
	fi
	log=$PWD/logs/$name.%A.log
	
	sbatch --partition=quick -D `pwd` --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=quick -D `pwd` --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
#fi

