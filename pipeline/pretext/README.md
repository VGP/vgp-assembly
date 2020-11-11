# Visualizing Hi-C with PretextMap and PretextView

Run Arima mapping pipeline + PretextMap to generate input file for PretextView.


## Requirements
* [PretextMap](https://github.com/wtsi-hpag/PretextMap) for converting `.bam` to `.pretext`
* [PretextView](https://github.com/wtsi-hpag/PretextView) for visualization, where you can visualize

## Workflow
```
Usage: ./_submit_pretext.sh <asm.fasta> <fastq.map> [jobid]
```

### Input
* `asm.fasta` : Reference to align. Symlink in the directory where this script is running.
* `fastq.map` : R1.fastq.gz <tab> R2.fastq.gz in one line. Concatenate before running.

### BWA Indexing
This step will be skipped if `$asm.fasta.fai` exists.

The first step generates bwa index of the `asm.fasta`.

```
index.sh $asm.fasta
```

### Arima mapping pipeline
This step will be skipped if `$asm.bam` exists. 

This code is a slightly modified version of the [ArimaGenomics mapping_pipeline](https://github.com/ArimaGenomics/mapping_pipeline).
```
arima_mapping_pipeline.sh fastq.map $asm $asm.fasta tmp [bwa_opts]
```

I have commented the last dedup step. Uncomment if you prefer to do the deduplication, run with picard tools.

### PretextMap
Convert `$asm.bam` to `$asm.pretext`

```
pretext.sh $asm.bam
```

### Output
* `$asm.pretext` : Input file for PretextView
