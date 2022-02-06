#!/bin/bash

#SBATCH --nodes=1
#SBATCH --ntasks=1

SAMPLE_LIST=($(<$1))
SAMPLE=${SAMPLE_LIST[${SLURM_ARRAY_TASK_ID}-1]}

echo "processing $SAMPLE"
echo "tmp dir ${TMPDIR}"

if [[ ${4} == "10x" ]] && [[ "$SAMPLE" == *"R1"* ]]; then

arg="-bc23"

fi

echo "\
time FastK -v -T${3} -N${2}/$(basename ${SAMPLE} .fastq.gz) -k${5} -t1 $SAMPLE -P${TMPDIR} $arg"
time FastK -v -T${3} -N${2}/$(basename ${SAMPLE} .fastq.gz) -k${5} -t1 $SAMPLE -P${TMPDIR} $arg
