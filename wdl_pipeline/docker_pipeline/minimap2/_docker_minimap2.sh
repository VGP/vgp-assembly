#!/bin/bash

if [[ -z $1 ]] || [[ ! -e input.fofn ]] ; then
        echo "Usage: ./_submit_minimap2.sh <ref.fasta>"
        echo -e "\tAssumes <input.fofn> in the same path"
	exit -1
fi

ref=$1
ref=`echo $ref | sed 's/.fasta$//g' | sed 's/.fa$//g' | sed 's/.fasta.gz$//g' | sed 's/.fa.gz$//g' `
#wait_for=$2

mkdir -p logs

# variable overwriting
export VGP_PIPELINE=/root/scripts

# index ref (if not present)
if [ ! -e $ref.idx ]; then
	cpus=16
	mem=32g
	name=index.$ref
	script=$VGP_PIPELINE/minimap2/minimap2_idx.sh
	args=$1
	log=logs/$name.log
	
	export SLURM_CPUS_PER_TASK=16
	echo -e "\n\n/bin/bash $script $args 2>&1 >$log"
	/bin/bash $script $args 2>&1 >$log

	#echo "\
	#sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args"
	#sbatch --partition=quick -D `pwd` --job-name=$name --time=240 --cpus-per-task=$cpus --mem=$mem --error=$log --output=$log $script $args > idx_jid
	#wait_for="--dependency=afterok:"`cat idx_jid`
fi

# align
if [ ! -e $ref.bam ]; then
	LEN=`wc -l input.fofn | awk '{print $1}'`
	cpus=16
	mem=32g
	name=minimap2.$ref
	script=$VGP_PIPELINE/minimap2/minimap2.sh
	args=$1

	export SLURM_CPUS_PER_TASK=16
	for i in $(seq 1 $END); do 
		log=logs/$name.$i.log
		echo -e "\n\n/bin/bash $script $args $i 2>&1 >$log"
		/bin/bash $script $args $i 2>&1 >$log
	done

	#echo "\
	#sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args"
	#sbatch --job-name=$name -a 1-$LEN $wait_for --cpus-per-task=$cpus --mem=$mem --partition=norm -D `pwd` --time=1-0 --error=$log --output=$log $script $args > minimap2_jid
	
	cpus=32
	mem=16g
	name=mrg.$ref
	script=$VGP_PIPELINE/minimap2/merge.sh
	args=$ref
	#walltime=12:00:00
	#dependency=""
	#if [ -e minimap2_jid ]; then
	#        dependency=`cat minimap2_jid`
	#        dependency="--dependency=afterok:$dependency"
	#fi
	log=$PWD/logs/$name.log

	echo -e "\n/bin/bash $script $args 2>&1 >$log"
        /bin/bash $script $args 2>&1 >$log

	#echo "\
	#sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	#sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
fi

