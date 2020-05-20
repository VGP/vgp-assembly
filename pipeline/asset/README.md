# Asset

## Installation

Download [asset](https://github.com/dfguan/asset)

Set `$asset` path variable pointing to asset top level directory
```
export asset=/path/to/asset
```

Scripts here are providing wrappers to run asset on the following platforms using input files:
* PacBio reads: runs minimap2 indexing and alignment
* 10XG Linked reads: runs bwa mem indexing and alignment, alternatively accepts bam file generated with longranger
* Illumina WGS reads: bwa mem indexing and alignment
* Bionano .bnx and .cmaps: Bionano Solve RefAligner (currently using Solve 3.3)
* Hi-C reads: runs bwa mem indexing and alignment, alternatively accepts bam file (ex. from Arima mapping pipeline)

`_submit*.sh` are slurm specific submitter scripts.
Other `*.sh` scripts can be run independently.

