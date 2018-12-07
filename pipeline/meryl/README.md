# K-mer counts with meryl

## Installation
Get [canu1.7](https://github.com/marbl/canu/releases/tag/v1.7) or [canu1.7.1]( https://github.com/marbl/canu/releases/tag/v1.7.1) release

*NOTE*: meryl in latest canu (1.8 or higher) is not compatible with this script. This script will be updated in the future to use the latest meryl. 

# Workflow
Use k=31 for obtaining the estimated genome size.
`_submit_10x_trimbc_meryl.sh` will submit jobs with its dependency, assuming we are running on slurm.

## Trim barcodes with scaff10x
On each 10X read fastq file, the barcodes will be trimmed and the headers will be renamed.

```
scaff_BC-reads-1 R1.fastq $genome_id.r1.$i.fastq read.name
scaff_BC-reads-2 $out_prefix.name R2.fastq $genome_id.r2.$i.fastq
```

## Concatenate re-headered fastq files
These files will be used for scaff10x, so let’s save them somewhere for later use.
```
cat $genome_id.r1.*.fastq | gzip > $genome_id.read-BC1.fastq.gz
cat $genome_id.r2.*.fastq | gzip > $genome_id.read-BC2.fastq.gz
```

## Build meryl db
After converting each fastq file to a simple fasta file, meryl starts creating k-mer db.
```
zcat $fastq | awk -v name=$name 'BEGIN {print ">"name".0"; num=1} {if (NR%4==2) print $1"N"; if (NR%4000000==0) {print "\n>"name"."num; num++}}' > $name.fa
meryl -B -C -s $name.fa -m $k -threads $SLURM_CPUS_PER_TASK -segments $SLURM_CPUS_PER_TASK -o $name.31
```
At the end, `$name.31.mcdat` and $name.31.mcidx` will be created.

## Merge meryl db and get histogram, filtered set, and simple stats
Merge all the meryl db sets:
```
meryl -M add -s $name1 -s $name2 … -o $genome_id.k31
```

Histogram:
```
meryl -Dh $genome_id.k31 > $genome_id.k31.hist
```

Filtering:
```
java -jar -Xmx512m kmerHistToPloidyDepth.jar $genome_id.k31.hist > $genome_id.k31.hist.ploidy
x=`sed -n 2p $genome_id.k31.hist.ploidy | awk '{print $NF}'`
meryl -M greaterthan $x -s $genome_id.k31-o $genome_id.k31.filt
```
Here we get the histogram and find the first lowest point.
The `$genome.k31.hist.ploidy` file is handy for getting the exact counts of haploid / diploid peaks and boundaries.

Estimated euchromatic region:
```
meryl -Dc -s $genome_id.k31.filt
```
The \‘Distinct\’ mers are k-mers appearing at least x times (filtering value $x above).
Now we are ready to run GenomeScope on our `$genome_id.k31.hist`.

# GenomeScope
The hist contains 4 columns, which is not accepted in the original GenomeScope.
For convenience, here is a [forked version of GenomeScope](https://github.com/arangrhie/genomescope/blob/master/genomescope.R) that runs locally.

Example:
```
Rscript genomescope.R $genome_id.k31.hist 31 120 $genome_id $genome_id 1
```

Will generate a foler named $genome_id, and plot the results with the title `$genome_id`.

## What to do next?
1.	Share the estimated genome size and heterozygosity with the vgp-assembly group.
2.	Upload `$genome_id.k31.*` to `s3://genomeark/species/<species_name>/<genome_id>/assembly_vgp_standard_1.5/intermediates/meryl/`.
3.	Upload `$genome_id.read-BC*.fastq.gz` to `s3://genomeark/species/<species_name>/<genome_id>/assembly_vgp_standard_1.5/intermediates/scaff10x/`.
4.	Perhaps, remove the intermediate files to save some disc-space.



