#!/bin/sh

#  Discover the job ID to run, from either a grid environment variable and a
#  command line offset, or directly from the command line.
#
if [ x$SLURM_ARRAY_TASK_ID = x -o x$SLURM_ARRAY_TASK_ID = xundefined -o x$SLURM_ARRAY_TASK_ID = x0 ]; then
  baseid=$1
  offset=0
else
  baseid=$SLURM_ARRAY_TASK_ID
  offset=$1
fi
if [ x$offset = x ]; then
  offset=0
fi
if [ x$baseid = x ]; then
  echo Error: I need SLURM_ARRAY_TASK_ID set, or a job index on the command line.
  exit
fi
jobid=`expr $baseid + $offset`
if [ x$SLURM_ARRAY_TASK_ID = x ]; then
  echo Running job $jobid based on command line options.
else
  echo Running job $jobid based on SLURM_ARRAY_TASK_ID=$SLURM_ARRAY_TASK_ID and offset=$offset.
fi

NUM=`wc -l input.fofn |awk '{print $1}'`
MAP_PRIMARY=0

if [ $jobid -gt $NUM ]; then
   echo "Error: id $jobid is out of bounds, max is $NUM"
   exit
fi

input=`head -n $jobid input.fofn |tail -n 1`
prefix=`echo $input |awk -F "/" '{print $NF}' |sed s/.fastq.gz//g |sed s/.fastq//g | sed s/.fasta.gz//g | sed s/.fasta//g`

ref="asm.fasta"
echo "Running with $input $prefix $ref"

module load samtools

if [ ! -e $prefix.sorted.bam ]; then
   if [ ! -e $prefix.sam ]; then
      echo "$tools/ngmlr/ngmlr-0.2.6/ngmlr -r $ref -q $input -t $SLURM_CPUS_PER_TASK -x pacbio --skip-write > $prefix.sam"
      $tools/ngmlr/ngmlr-0.2.6/ngmlr -r $ref -q $input -t $SLURM_CPUS_PER_TASK -x pacbio --skip-write > $prefix.sam
   fi
   echo "samtools sort -@$SLURM_CPUS_PER_TASK -O bam -o $prefix.sorted.bam -T $prefix.tmp $prefix.sam" &&
   samtools sort -@$SLURM_CPUS_PER_TASK -O bam -o $prefix.sorted.bam -T $prefix.tmp $prefix.sam &&
   echo "samtools index $prefix.sorted.bam" &&
   samtools index $prefix.sorted.bam &&
   echo "rm $prefix.sam" &&
   rm $prefix.sam ||
   echo "[ ERROR ] :: encounterd an issue while running $prefix.sam"
else
   echo "Already done!"
fi
