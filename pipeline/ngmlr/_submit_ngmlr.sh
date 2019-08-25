#! /bin/bash

if [ ! -e "input.fofn" ] || [ ! -e "asm.fasta" ]; then
	echo "Run ngmlr in parallel. Requires <input.fofn> and <asm.fasta>"
	echo "No input.fofn or asm.fasta found. Exiting."
	exit -1;
fi

if [ -z $1 ]; then
	echo "Usage: ./run.sh <jobname>"
	exit -1;
fi

cpus=32
mem=52g
name=ngmlr_$1
script=$VGP_PIPELINE/ngmlr/ngmlr.sh
walltime=3-0

mkdir -p logs
log=$PWD/logs/$name.%A_%a.log

LEN=`wc -l input.fofn | awk '{print $1}'`

if ! [ -e $1.bam ]; then
	arr=""
	for i in $(seq 1 $LEN)
	do
		prefix=`sed -n ${i}p input.fofn`
		prefix=`echo $prefix |awk -F "/" '{print $NF}' |sed s/.fastq.gz//g |sed s/.fastq//g | sed s/.fasta.gz//g | sed s/.fasta//g `
		if ! [ -e $prefix.sorted.bam ]; then
			arr=$arr${i},
		fi
	done
	echo $arr
	arr_count=`echo $arr | tr -cd , | wc -c`

	if [ "$arr_count" -eq 0 ]; then
		echo "sorting finished. skip alignment."
		rm -f $name.jid
	elif [ "$arr_count" -lt 20 ]; then
		echo $arr_count
		echo "sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --array=$arr --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
		sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --array=$arr --mem=$mem --time=$walltime --error=$log --output=$log $script $args > $name.jid

	else
		echo "sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --array=1-$LEN --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
		sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --array=1-$LEN --mem=$mem --time=$walltime --error=$log --output=$log $script $args > $name.jid
	fi
fi

cpus=32
mem=16g
name=mrg_$1
script=$VGP_PIPELINE/ngmlr/merge.sh
args=$1
walltime=2-0
dependency=""
if [ -e ngmlr_$1.jid ]; then
	dependency=`cat ngmlr_$1.jid`
	dependency="--dependency=afterok:$dependency"
fi
log=$PWD/logs/$name.%A.log
echo "\
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name $dependency --mem=$mem --time=$walltime --error=$log --output=$log $script $args
