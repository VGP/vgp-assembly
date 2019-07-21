#!/bin/bash

if [[ -z $1 ]] || [[ ! -e input.fofn ]] ; then
	echo "Usage: ./_submit_purge_dups.sh <asm.fasta> <input.fofn> [jid]"
	echo -e "\t[jid]: set --dependency=afterok:[jid]"
	exit -1
fi

asm=$1
ref=`echo $asm | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' `

fofn=$2
wait_for=$3
if [ ! -z $wait_for ]; then
	wait_for="--dependency=afterok:"$wait_for
fi

mkdir -p logs

if [ ! -e $ref.idx ]; then
	cpus=16
	mem=32g
	name=minimap2.index.$ref
	script=$VGP_PIPELINE/purge_dups/minimap2_idx.sh
	args=$1
	log=logs/$name.%A.log
	echo "\
	sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
	sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args > idx_jid
	wait_for="--dependency=afterok:"`cat idx_jid`
fi

if [ ! -s purged.fa ]; then
	LEN=`wc -l $fofn | awk '{print $1}'`
	NUM_PAF=`ls *.read.paf.gz 2> /dev/null | wc -l | awk '{print $1}'`
	echo "$NUM_PAF pafs found."
	if [[ $LEN -gt $PAFs ]]; then
		cpus=16
		mem=32g
		name=minimap2.map.$ref
		script=$VGP_PIPELINE/purge_dups/minimap2.sh
		args=$1
		log=logs/$name.%A_%a.log
	
		echo "\
		sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args"
		sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args > minimap2_jid
	fi

        if [[ ! -s $asm.split.self.paf.gz ]]; then
                cpus=16
                mem=32g
                name=minimap2_self.$ref
                script=$VGP_PIPELINE/purge_dups/minimap2_self.sh
                args=$asm
                log=logs/$name.%A.log
		wait_for=""

                echo "\
                sbatch --job-name=$name $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args"
                sbatch --job-name=$name $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args > self_jid
        fi

	cpus=4
	mem=12g
	name=purge_dups.$ref
	script=$VGP_PIPELINE/purge_dups/purge_dups.sh
	args=$asm
	walltime=1:00:00
	dependency=""
	if [[ -e minimap2_jid || -e self_jid ]]; then
	        dependency=`cat minimap2_jid`
		dependency=$dependency,`cat self_jid`
	        dependency="--dependency=afterok:$dependency"
	fi

	log=$PWD/logs/$name.%A.log
	echo "\
	sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
fi

