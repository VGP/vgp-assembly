#!/bin/bash

module load bwa/0.7.17
module load gcc/7.4.0

name=$1

# -nodes: Reserve 2 cores for computational overheads

echo "
$tools/scaff10x/Scaff10X-4.1/src/scaff10x -nodes $((SLURM_CPUS_PER_TASK-2)) -longread 1 -gap 100 -matrix 2000 -read-s1 12 -read-s2 8 -link 10 -score 20 -edge 50000 -block 50000 -data input.dat asm.fasta $name.scaff10x.fasta
"
$tools/scaff10x/Scaff10X-4.1/src/scaff10x -nodes $((SLURM_CPUS_PER_TASK-2)) -longread 1 -gap 100 -matrix 2000 -read-s1 12 -read-s2 8 -link 10 -score 20 -edge 50000 -block 50000 -data input.dat asm.fasta $name.scaff10x.fasta

