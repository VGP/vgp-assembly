# Freebayes
For installing and running `longranger`, see [this page](https://github.com/VGP/vgp-assembly/blob/master/pipeline/longranger).

## Requirements
* Longranger: see [this page](https://github.com/VGP/vgp-assembly/blob/master/pipeline/longranger) for installation and setting up input fastq path
* [FreeBayes v1.3.1](https://github.com/ekg/freebayes/releases)
* Bcftools v1.8+
* Samtools


## Workflow
This folder contains workflows for
* `_submit_longranger_freebayes.sh` : Longranger + FreeBayes + Consensus + QV
* `_submit_longranger_freebayes_large.sh` : Longranger + FreeBayes + Consensus + QV for large genomes (Recommended for genomes > 4Gbp)
* `_submit_freebayes.sh` : FreeBayes + Consensus + QV, for running directly on a bam file
* `_submit_bwa_freebayes.sh`: BWA + FreeBayes + Consensus + QV, for regular Illumina WGS data

### FreeBayes
Assuming `longranger align` finished with no issues, we will have `$sample/outs/possorted_bam.bam` where `$sample` is the longranger output directory.

The submitter scripts links the following longranger output files:
```
	ln -s $sample/outs/possorted_bam.bam aligned.bam
	ln -s $sample/outs/possorted_bam.bam.bai aligned.bam.bai
	ln -s $sample/outs/summary.csv
```

And assumes that longranger generated `genome.fa.fai` under `refdata-asm/fasta/`.

```
freebayes_v1.3.1.sh $sample [array_idx]
```
This script is for running freebayes in parallel, for every other 100th scaffold in 100 array jobs.
The following can be run on a single instance with multiple threads, or on multiple machines:

```shell
for array_idx in $(seq 1 100)
do
  ./freebayes_v1.3.1.sh $sample $array_idx
done
```

Mean coverage (mean_cov) is extracted from `summary.csv`. If not, default value of 50 will be used.
Regions with >600x (50 x 12) depth of coverage will be ignored.

The actual code for running freebayes is as following:
```
freebayes --bam $bam --region=$contig:1-$end --skip-coverage $((mean_cov*12)) -f $fasta | bcftools view --no-version -Ou > bcf/$contig_no_pipe.bcf
```

### Consensus
Concatenate and normalise non-REF variants
```
bcftools concat -f concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o $sample.bcf --threads $threads
```

Filter: This will do very basic filtering of the callset (`QUAL>1`) and select only homozygous ALT and heterozygous non-REF sites.
```
bcftools view -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Oz --threads=$threads $sample.bcf > $sample.changes.vcf.gz
bcftools index $sample.changes.vcf.gz
```

Consensus: At heterozygous non-REF sites (both alleles do not match the reference fasta), the longest allele will be chosen.
```
bcftools consensus -Hla -f $fasta $sample.changes.vcf.gz > ${sample}_fb.fasta
```

Alternatively, if you have a phased VCF and would like to consistently pick one of the haplotypes at all heterozygous sites, set `-H` to `1` or `2`:
```
bcftools consensus -H1 -f $fasta $sample.changes.vcf.gz > $sample.fasta
```

Get the number of bases subject to be polished
```
bcftools view -H -i 'QUAL>1 && (GT="AA" || GT="Aa")' -Ov $sample.changes.vcf.gz | awk -F "\t" '{print $4"\t"$5}' | awk '{lenA=length($1); lenB=length($2); if (lenA < lenB ) {sum+=lenB-lenA} else if ( lenA > lenB ) { sum+=lenA-lenB } else {sum+=lenA}} END {print sum}' > $sample.numvar
echo "Num. bases affected: `cat $sample.numvar`"
```

### Mapping based QV

Genome coverage
```
samtools view -F 0x100 -u $bam | bedtools genomecov -ibam - -split > aligned.genomecov
```

Bases excluding regions with too low (<3x) or too high (> mean_cov x 12) coverage
```
awk -v l=$l_filter -v h=$h_filter '{if ($1=="genome" && $2>l && $2<h) {numbp += $3}} END {print numbp}' aligned.genomecov > $genome.numbp
```

Calculate QV
```
NUM_VAR=`cat $sample.numvar` # From the above
echo "Total num. bases subject to change: $NUM_VAR"
QV=`echo "$NUM_VAR $NUM_BP" | awk '{print (-10*log($1/$2)/log(10))}'`
echo $QV > $sample.qv
echo "QV of this genome $sample: $QV"
```
The actual script for calculating QV is in a [different directory](https://github.com/VGP/vgp-assembly/blob/master/pipeline/qv) and can be run separately. Note the mapping based QV is a very rough measure, and is sensitive to coverage fluctuation and variant calling as discussed in our [VGP paper](https://doi.org/10.1101/2020.05.22.110833).

### Troubleshooting
* Failed freebayes job: 
For a failed job, remove the failed scaffold’s .bcf under `bcf` dir and restart `freebayes_v1.3.sh` or launch `_submit_freebayes.sh`.
It will skip if bcf/scaffold.bcf is already present in the directory and will only run for those that aren’t generated.
