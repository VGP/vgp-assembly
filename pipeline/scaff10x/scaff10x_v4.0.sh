#!/bin/bash

module load bwa/0.7.17
module load gcc/7.4.0

name=$1

echo "
/data/Phillippy/tools/scaff10x/Scaff10X-4.0/src/scaff10x -nodes $SLURM_CPUS_PER_TASK -longread 1 -gap 100 -matrix 2000 -read-s1 12 -read-s2 8 -link 10 -score 20 -edge 50000 -link 8 -block 50000 -data input.dat asm.fasta $name.scaff10x.fasta
"
/data/Phillippy/tools/scaff10x/Scaff10X-4.0/src/scaff10x -nodes $SLURM_CPUS_PER_TASK -longread 1 -gap 100 -matrix 2000 -read-s1 12 -read-s2 8 -link 10 -score 20 -edge 50000 -link 8 -block 50000 -data input.dat asm.fasta $name.scaff10x.fasta

