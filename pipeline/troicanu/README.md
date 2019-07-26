# Trio canu

Obtain the v1.9 branch

```
git clone --single-branch --branch v1.9 https://github.com/marbl/canu.git
```

## Install
```
cd canu/src
make -j <number of threads>
```

## Run

Here is an example of how to run a 3.1G genome.
This script assumes we have `paternal` and `maternal` dir with illumina reads and `pacbio` dir for the PacBio subreads of the child.

*NOTE*: VGP subreads comes in bam file. A simple (so far the fastest) way to convert a bam file to fasta file is:
```
samtools view subreads.bam | awk '{print ">"$1"\n"$10}' | pigz -c > subreads.fasta.gz
```

Following is the trio canu options run for the `assembly_nhgri_trio_1.6` pipeline:

```
export PATH=/path/to/canu/bin:$PATH

canu -p asm -d triocanu genomeSize=3.1g \
    -haplotypePaternal paternal/*.fastq.gz \
    -haplotypeMaternal maternal/*.fastq.gz \
    -pacbio-raw pacbio/*subreads.fasta.gz \
    correctedErrorRate=0.105 \
    hapUnknownFraction=0.01 \
    corMhapSensitivity=normal \
    corMinCoverage=0 \
    useGrid=true gridOptions="--time=3-0" \
    hapMemory=48 batMemory=128g \
    stageDirectory="/lscratch/\$SLURM_JOBID" \
    gridEngineStageOption="--gres=lscratch:100"
```

* `correctedErrorRate=0.105`: Tuned for Sequel1 instruments, where most of our VGP genomes are sequenced with
* `hapUnknownFraction=0.01`: Include unknown reads (not enough hap-specific markers supported to be binned) for assembly when unknown base fraction is > 0.01.
* `corMhapSensitivity=normal` : For genomes sequenced with ~60x, binned reads marginally end up to be less than 30x. By default, MHAP is set to `sensitive` for low coverage sequenced genomes. To avoid sequence noise being targeted for overlaps, which causes frequent breaks in the assembly, set this slightly less sensitive.
* `corMinCoverage=0` : Use all reads.
* `hapMemory` and `batMemory` needs to be set.

The rest are tuned for the slurm environment local at NIH. `lscratch` is the temporary local scratch disk to prevent overhead IO for shared environments.

# Other options

For Sequel 2 data, use the following options (Courtesy by Sergey Koren):
`correctedErrorRate=0.035 utgOvlErrorRate=0.065 trimReadsCoverage=2 trimReadsOverlap=500`

For Nanopore guppy flip-flop HAC:
`correctedErrorRate=0.105 corMhapOptions="--threshold 0.8 --num-hashes 512 --ordered-sketch-size 1000 --ordered-kmer-size 14"`

# Useful outputs

Besides the contigs.fasta and .report files, following files can tell the haplotype binning performance and the expected sequence coverage used while canu is running.

* `triocanu/haplotype/haplotype.log` : Contains the num. of reads and bases in each bins.
* `triocanu/haplotype/0-meryl/*.meryl` : Can be used for evaluating hap-specific kmers in the final assembly set.
Needs to be filtered for low coverage erroneous k-mers.


