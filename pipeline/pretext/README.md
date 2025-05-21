# Visualizing Hi-C with PretextMap and PretextView

Run Arima mapping pipeline + PretextMap to generate input file for PretextView.


## Requirements
* [PretextMap](https://github.com/wtsi-hpag/PretextMap) for converting `.bam` to `.pretext`
* [PretextView](https://github.com/wtsi-hpag/PretextView) for visualization, where you can visualize

## Workflow
```
Usage: ./_submit_pretext.sh <asm.fasta> <fastq.map> <out_prefix> [jobid]
```

### Input
* `asm.fasta`  : Reference to align. Symlink in the directory where this script is running.
* `fastq.map`  : R1.fastq.gz <tab> R2.fastq.gz in one line. Concatenate before running.
* `out_prefix` : Output prefix. Final output will be a merged <out_prefix>.bam and <out_prefix>.pretext

### BWA Indexing
This step will be skipped if `$asm.fasta.fai` exists.

The first step generates bwa index of the `asm.fasta`.

```
index.sh $asm.fasta
```

### Arima mapping pipeline
This step will be skipped if `$asm.bam` exists. 

This code is a slightly modified version of the [ArimaGenomics mapping_pipeline](https://github.com/ArimaGenomics/mapping_pipeline).
`arima_mapping_pipeline.sh` is under [salsa](../salsa/arima_mapping_pipeline.sh).
```sh
arima_mapping_pipeline.sh fastq.map $out_prefix $asm.fasta tmp [bwa_opts]
```

I have commented the last dedup step. Uncomment if you prefer to do the deduplication, run with picard tools.

### PretextMap
Convert `$out_prefix.bam` to `$out_prefix.pretext`

```sh
pretext.sh $out_prefix.bam
```

### Output
* `$out_prefix.bam`     : In case further inspection is needed. Note it's sorted by read ID.
* `$out_prefix.pretext` : Input file for PretextView


### Update
Confirmed `.pretext` is compatable with [PretextViewAI v1.0.3 release](https://github.com/sanger-tol/PretextView/releases/tag/1.0.3).
