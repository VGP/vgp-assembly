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

NUM_JOBS=`wc -l input.fofn |awk '{print $1}'`

if [ $jobid -gt $NUM_JOBS ]; then
  echo Error: Only $NUM_JOBS partitions, you asked for $jobid.
  exit 1
fi

jobid=`echo $jobid |awk '{print $1}'`
input_file=`head -n $jobid input.fofn |tail -n 1`

if [[ "$input_file" =~ \.bam$ ]] ; then
    if [[ ! -e "${input_file/.bam/.fasta}" ]] ; then
        module load samtools
       echo "convert bam to fasta"
        samtools view $input_file | awk '{print ">"$1"\n"$10}' > ${input_file/.bam/.fasta}
    fi
    input_file=${input_file/.bam/.fasta}
fi

output_file=`echo $input_file |awk -F "/" '{print $(NF)}' |sed s/.fasta.gz/.out/g |sed s/.fastq.gz/.out/g |sed s/.fasta/.out/g |sed s/.fastq/.out/g`

echo $input_file
echo $output_file

if [ ! -s $output_file ]; then
   module load mash
   mash screen -w $VGP_PIPELINE/mash/refseq.genomes.k21s1000.msh $input_file > $output_file
fi
