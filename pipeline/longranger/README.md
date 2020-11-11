# Longranger
For polishing with longranger and freebayes, go to [this page](https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish) and launch `_submit_longranger_freebayes.sh`.


This directory contains individual scripts to run `longranger align`.

## Requirements
* [Longranger-2.2.2](https://support.10xgenomics.com/genome-exome/software/pipelines/latest/installation) installed as `$tools/longranger/longranger-2.2.2/longranger`

Modify template accordingly, under `longranger-2.2.2/martian-cs/2.3.2/jobmanagers` for lsf, pbspro, sge, slurm, or torque.

Below is how my `slurm.template` looks like:
```#!/usr/bin/env bash
#
# Copyright (c) 2016 10x Genomics, Inc. All rights reserved.
#
# =============================================================================
# Setup Instructions
# =============================================================================
#
# 1. Add any other necessary Slurm arguments such as partition (-p) or account
#    (-A). If your system requires a walltime (-t), 24 hours (24:00:00) is
#    sufficient.  We recommend you do not remove any arguments below or Martian
#    may not run properly.
#
# 2. Change filename of slurm.template.example to slurm.template.
#
# =============================================================================
# Template
# =============================================================================
#
#SBATCH -J __MRO_JOB_NAME__
#SBATCH --export=ALL
#SBATCH --cpus-per-task=__MRO_THREADS__
#SBATCH --signal=2 -p norm -t 1-0
#SBATCH --no-requeue
### Alternatively: --ntasks=1 --cpus-per-task=__MRO_THREADS__
###   Consult with your cluster administrators to find the combination that
###   works best for single-node, multi-threaded applications on your system.
#SBATCH --mem=__MRO_MEM_GB__G
#SBATCH -o __MRO_STDOUT__
#SBATCH -e __MRO_STDERR__

__MRO_CMD__
```
The most important part is the `-p norm -t 1-0` part, which will be depending on your job scheduler configuration.

## Workflow
```
_submit_longranger.sh <genome> <ref.fasta>
```

### Index
```
$tools/longranger/longranger-2.2.2/longranger mkref $ref.fasta
```

### Align
`longranger.sh` runs the following:
```
$tools/longranger/longranger-2.2.2/longranger align \
--id=$genome \
--fastq=/path/to/genomic_data/10x/ \
--sample=$genome \
--reference=refdata-$ref \
--jobmode=slurm \
--localcores=32 \
--localmem=60 \
--maxjobs=500 \
--jobinterval=5000 \
--disable-ui \
--nopreflight
```

Note the path to `--fastq` is fixed. Replace that to your file path, accordingly.
The `10x/` directory expects to have the input fastq files in this format:
* <prefix>_S*_L00*_R1_001.fastq.gz
* <prefix>_S*_L00*_R2_001.fastq.gz

For large genomes, where the reference asm.fasta is larger than 4 Gbp, use the `longranger_4G.sh` instead.
This uses an [override](https://github.com/VGP/vgp-assembly/blob/master/pipeline/longranger/override_4G.json) configuration to martian-cs, with larger memory.

### Output
* `$genome/outs/possorted_bam.bam` and `.bai`
* `refdata-asm/fasta/genome.fa.fai` and other indexes
