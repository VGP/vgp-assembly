#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./meryl_build.sh <k-size> <input.fofn> [num_line]"
	exit -1
fi

num_line=$SLURM_ARRAY_TASK_ID 
if [ -z "$num_line" ] ; then
	num_line=$3
fi

k=$1
fastq=`sed -n ${num_line}p $2 | awk '{print $1}'`
name=`echo $fastq | sed 's/.fasta.gz$//g' | sed 's/.fq.gz$//g' | sed 's/.fasta$//g' | sed 's/.fa$//g'`

if [ -e "$name.$k.mcdat" ]; then
	echo "$name.$k.mcdat already exists. skip building."
	exit 0;
fi

echo "PWD : $PWD"
if ! [[ -e $name.fa ]] || ! [[ -e $name.fasta ]]; then
	echo "convert $fastq to FASTA"
	if [[ ${fastq} =~ \.gz$ ]]; then
		zcat $fastq | awk -v name=$name 'BEGIN {print ">"name".0"; num=1} {if (NR%4==2) print $1"N"; if (NR%4000000==0) {print "\n>"name"."num; num++}}' > $name.fa
	else
		cat $fastq | awk -v name=$name 'BEGIN {print ">"name".0"; num=1} {if (NR%4==2) print $1"N"; if (NR%4000000==0) {print "\n>"name"."num; num++}}' > $name.fa
	fi
fi

echo "build meryl db"
echo "\
$tools/canu/canu-1.7/Linux-amd64/bin/meryl -B -C -s $name.fa -m $k -threads $SLURM_CPUS_PER_TASK -segments $SLURM_CPUS_PER_TASK -o $name.$k"
$tools/canu/canu-1.7/Linux-amd64/bin/meryl -B -C -s $name.fa -m $k -threads $SLURM_CPUS_PER_TASK -segments $SLURM_CPUS_PER_TASK -o $name.$k
#rm $name.fa
#rm $name.fa.fastaidx
