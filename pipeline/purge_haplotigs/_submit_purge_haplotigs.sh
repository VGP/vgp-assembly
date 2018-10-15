#! /bin/bash

if [ -z $1 ]; then
	echo "Usage: ./_submit_purge_haplotigs.sh <name>"
	echo -e "\tRequires aligned.bam"
	exit -1
fi

mkdir -p logs

if ! [ -e aligned.bam ]; then
	echo "Requires aligned.bam. Exit."
	exit -1
fi

waitfor=""
if ! [ -e coverage_stats.csv ]; then
	cpus=4
	mem=8g
	name=${1}_step1_2
	script=$VGP_PIPELINE/purge_haplotigs/_step1_2.sh
	args=$1
	walltime=2-0

	log=logs/$name.%A_%a.log

	echo "\
	sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > $name.jid
	waitfor=`cat $name.jid | tail -n1`
	waitfor="--dependency=afterok:$waitfor"
fi

if ! [ -e ${1}_curated.fasta ]; then
        cpus=32
        mem=80g
        name=${1}_step3
	partition=norm
        script=$VGP_PIPELINE/purge_haplotigs/_step3.sh
        args=$1
        walltime=5-0

        log=logs/$name.%A_%a.log

        echo "\
        sbatch --partition=$partition $waitfor --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
        sbatch --partition=$partition $waitfor --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args 
fi
