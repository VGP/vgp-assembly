# Assembly pipeline for the VGP genomes

This folder contains wrappers to run individual tools in each assembly / evaluation step.

## Pipeline

### Before running an assembly (QC)
* [Meryl k-mer counting for genome size, repeat contents, heterozygosity estimates](https://github.com/VGP/vgp-assembly/tree/master/pipeline/meryl)
* [Mash Distance and Screen](https://github.com/VGP/vgp-assembly/tree/master/pipeline/mash)

### Contigging
* FALCON-Unzip: FALCON-Unzip was run on DNAnexus. See [this page](https://github.com/VGP/vgp-assembly/tree/master/dx_workflows/vgp_falcon_and_unzip_assembly_workflow) for more details.
* [TrioCanu](https://github.com/VGP/vgp-assembly/tree/master/pipeline/triocanu)

### Purging
* [Minimap2 + purge_dups](https://github.com/VGP/vgp-assembly/tree/master/pipeline/purge_dups)
* [NGMLR](https://github.com/VGP/vgp-assembly/tree/master/pipeline/ngmlr) + [purge_haplotigs](https://github.com/VGP/vgp-assembly/tree/master/pipeline/purge_haplotigs) _retired_

### Scaffolding
* [Scaff10x](https://github.com/VGP/vgp-assembly/tree/master/pipeline/scaff10x)

### Polishing
* [Bionano](https://github.com/VGP/vgp-assembly/tree/master/pipeline/bionano): Hybrid Scaffolding and other scripts to run bionano software
* Hi-C
  * [Arima mapping + SALSA](https://github.com/VGP/vgp-assembly/tree/master/pipeline/salsa)
  * Visualizing Hi-C map with [PretextMap and PretextView](https://github.com/VGP/vgp-assembly/tree/master/pipeline/pretext)
  * Visualizing Hi-C map with [Juicer](https://github.com/VGP/vgp-assembly/tree/master/pipeline/pretext) _retired_

### Polishing
* [Arrow polishing on Grid](https://github.com/skoren/ArrowGrid): Note this is an external repository
* [Longranger](https://github.com/VGP/vgp-assembly/tree/master/pipeline/longranger) + [Freebayes](https://github.com/VGP/vgp-assembly/tree/master/pipeline/freebayes-polish)
* [Bwa mem + FreeBayes](https://github.com/VGP/vgp-assembly/tree/master/pipeline/bwa): for polishing with no 10XG data
* [Mapping-based QV](https://github.com/VGP/vgp-assembly/tree/master/pipeline/qv): Get mappable regions and compute QV based on variant calls

### Post evaluation
* [Stats](https://github.com/VGP/vgp-assembly/tree/master/pipeline/stats): Basic assembly continuity - N-xx or NG-xx
* [Asset](https://github.com/VGP/vgp-assembly/tree/master/pipeline/asset): Reliable blocks
* [Telomere](https://github.com/VGP/vgp-assembly/tree/master/pipeline/telomere): Find telomere motives
* [Merqury](https://github.com/marbl/merqury): _k-mer_ based completeness, base pair QV, phase block statistics. Note this is an external repository
* [KAT](https://github.com/VGP/vgp-assembly/tree/master/pipeline/kat): spectra-cn with `kat comp` _retired_
* [BUSCO](https://github.com/VGP/vgp-assembly/blob/master/pipeline/busco): _warning_ this had multithreading issues and did not properly run
* [RepeatMasker](https://github.com/VGP/vgp-assembly/tree/master/pipeline/repeatmasker): _warning_ Not fully tested
* [MashMap2 alignment](https://github.com/VGP/vgp-assembly/tree/master/pipeline/mashmap): For quick and approximate genome-to-genome alignment
* [MUMMER nucmer alignment](https://github.com/VGP/vgp-assembly/tree/master/pipeline/nucmer): For exact genome-to-genome alignment

## Environments setting
Each script was written and tested to run on a slurm scheduler.
To begin, two environment variables are assumed:
* $VGP_PIPELINE=/path/to/this/git/vgp-assembly/pipeline
* $tools=/path/where/tools/are/installed

Add the path to the pipeline in `~/.bash_profile`
```
export VGP_PIPELINE=/path/to/git/vgp-assembly/pipeline
export tools=/path/where/tools/are/installed
```
And `source ~/.bash_profile.

Try `./_submit_` in your desired step to see required input files. Most of the scripts only require `.fasta` and `input.fofn` files as input.

The `_submit_` scripts are wrappers to submit jobs to launch the actual worker scripts.
In this file, required numbers of CPUs and memory usage in each step are specified as `cpus` and `mem`, optimized to run on 1 ~ 5 Gbp genomes. Some folders contain `_submit_*_large.sh`, which were optimized to run with genomes over 3 Gbp. Most jobs finish before the specified `walltime`, which is formatted in `hh:mm:ss` or `day-0`.

Most of the tools are loaded as modules. Adjust the module load part to your own environment.

List of tools assumed loadable or accessible with no path are:
* [Samtools](https://github.com/samtools/samtools)
* [Minimap2](https://github.com/lh3/minimap2)
* [BEDTools](https://github.com/arq5x/bedtools2)
* [BWA](https://github.com/lh3/bwa)

## VGP Assembly Paper
Input data and scripts for figures presented in the [VGP paper](https://doi.org/10.1101/2020.05.22.110833) are [in this repo](https://github.com/arangrhie/Scratch/tree/master/VGP).
