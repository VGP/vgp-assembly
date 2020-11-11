# Arima mapping + SALSA2

## Requirements
* [BWA](https://github.com/lh3/bwa)
* [Samtools](https://github.com/samtools/samtools)
* [BEDTools](https://github.com/arq5x/bedtools2)
* [Salsa2.2](https://github.com/machinegun/SALSA/releases/tag/v2.2) (See below)

### Salsa2 on a conda environment
Download [Salsa2.2](https://github.com/machinegun/SALSA/releases/tag/v2.2) in your $tools/salsa path

```
cd $tools
mkdir salsa2
wget https://github.com/machinegun/SALSA/archive/v2.2.tar.gz
tar -zxf v2.2.tar.gz
```

Build python environment: This needs to be done only once

```
mkdir -p $tools/conda/temp
cd $tools/conda/temp
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -p $tools/conda -b
source $tools/conda/etc/profile.d/conda.sh
conda update conda
conda create -n salsa_env python=2.7.15 networkx==1.11
```

Test if salsa is working
```
source $tools/conda/etc/profile.d/conda.sh
conda activate salsa_env
python $tools/salsa2/SALSA-2.2/run_pipeline.py -h
```

## Workflow
```
Usage: ./_submit_salsa_2.2.sh <fasta> [jobid]
```

### Input
* `fastq.map` : R1.fastq.gz <tab> R2.fastq.gz in one line. Concatenate before running.
* `asm.fasta` : Reference to align. Symlink in the directory where this script is running.
* `re_bases.txt` : Restriction enzyme site bases.

### BWA Indexing
The first step generates bwa index of the `asm.fasta`.
```
index.sh $asm.fasta
```

### Arima mapping pipeline
This code is a slightly modified version of the [ArimaGenomics mapping_pipeline](https://github.com/ArimaGenomics/mapping_pipeline).
```
arima_mapping_pipeline.sh fastq.map out_prefix $asm.fasta tmp [bwa_opts]
```

I have commented the last dedup step. Uncomment if you prefer to do the deduplication, run with picard tools.

### Salsa 2.2
Convert the bam to bed and finally do the scaffolding
```
salsa2.2.sh $asm out_prefix
```

### Output
* `out_prefix/SCAFFOLDED_FINAL.fasta`
