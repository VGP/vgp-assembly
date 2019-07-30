#!/bin/bash

set -e

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF

	bamextract.sh extracts reads from bam file(s) based on a text file with read names

	Required arguments are:
	-r the text file with the list of the reads to be extracted
	-b the folder where bams file from which reads need to be extracted can be found
	
	Optional arguments are:
	-p partition
	
	Picard is required. Please export the folder containing picard.jar file to the env variable PICARD_PATH

EOF

exit 0

fi

while getopts ":r:b:p:" opt; do

	case $opt in
		r)
			LIST=$OPTARG
			echo "Read names: -r $OPTARG"
			;;
		b)
			BAM=$OPTARG
			echo "Bam files folder: -b $OPTARG"
			;;
		p)
			PAR=$OPTARG
			echo "partition: -p $OPTARG"
			;;
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac

printf "\n"

done

printf "\n"

echo Starting at `date`
echo This is job $SLURM_JOB_ID
echo Running on `hostname`

# if ! [[ -e "reads.fasta" ]]; then
# 
# 	gunzip -c ${FASTA} > reads.fasta
# 
# fi

# if ! [[ -e "read_names.txt" ]]; then
# 
# 	grep -o -P '(?<=>)(\S*)' reads.fasta | sort | uniq -u > read_names.txt
# 
# fi

if ! [[ -e "bam.ls" ]]; then

	ls -1 ${BAM}/*.bam > bam.ls

fi

if ! [[ -e "read_names.uniq" ]]; then

	grep -o "m[0-9]*_[0-9]*_[0-9]*" ${LIST} | sort | uniq > read_names.uniq

fi

if ! [[ -e "bam.uniq" ]]; then

	grep -o "m[0-9]*_[0-9]*_[0-9]*" bam.ls | sort | uniq > bam.uniq

fi

chk1=`cksum read_names.uniq | awk -F" " '{print $1}'`
chk2=`cksum bam.uniq | awk -F" " '{print $1}'`

if [ $chk1 -eq $chk2 ]
then
  echo "All necessary bam files provided"
else
  echo "WARNING: one or more bam files appear to be missing!"
  exit 1
fi

printf "\n"

if ! [[ -e "reads/" ]]; then

	mkdir reads	

while read p; do

  echo "writing $p to list"
  
  grep ${p} ${LIST} | awk '{print $2}' > reads/$p.ls

done <read_names.uniq

fi

printf "\n"

if ! [[ -e "filtered/" ]]; then

	mkdir filtered	

fi

if ! [[ -e "logs/" ]]; then

	mkdir logs	

fi

if ! [ -z ${PAR} ]; then

	PAR="--partition=${PAR}"

fi

NFILES=$(cat bam.uniq | wc -l)	

sbatch ${PAR} --array=1-${NFILES} --cpus-per-task=1 -o logs/slurm-%A_%a.out extract.sh