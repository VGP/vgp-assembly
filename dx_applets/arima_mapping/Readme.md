# Arima mapping

## What does this app do?

This app mapps HiC paired-end reads obtained via Arima HiC kits or services to assembled contigs in order to scaffold the contigs. 

## What are typical use cases for this app?

This app can be used to map HiC reads from an Arima HiC kit or service to a set of assembled contigs, producing a bam file which 
can be fed to a downstream scaffolding algorithm to scaffold the contigs.

This app uses the BWA (Burrows-Wheeler Alignment Tool) software package, and in particular the BWA-MEM routine.  HiC reads, unlike 
traditional paired-end reads, can have wildly different genomic distances between the two paired reads.  This is due to the fact that
in HiC experiments, the sample is first crosslinked to preserve the 3D conformation of the genome.  This crosslinked DNA is then 
digested, and spatially proximal free ends of DNA are then ligated.  This process preserves both short-range and long-range information.
As a result, the frequency with which two pieces of DNA are seen involved in the same paired-end read can provide an indication of 
the genomic distance that separates the sequences.  

Because of the variable distance between the forward and reverse reads and because a given forward/reverse read can be a chimera of two 
genomic regions, special care must be used during mapping.  In particular, this routine will map each of the forward and reverse
reads separately using BWA-MEM.  These separately mapped files are then examined and potential chimeric reads are filtered and corrected
to retain only the 5' end of the chimera (since the 3' side can originate from the same contiguous DNA found as the 5' side of the pair.)

Once the two mapping files have been filtered to remove chimeric regions, they are combined, sorted, and duplicates are marked.  

## What data are required for this app to run?

This app requires reads files in gzipped FASTQ format (`*.fastq.gz` or `*.fq.gz`) generated from a HiC experiment using an Arima kit or service.

The app also requires a reference genome sequence index. This must be a gzipped tar archive file (`*.bwa-index.tar.gz`) containing
all the sequence index files as previously output by the BWA indexer. (Indexing is an one-time operation that needs to be performed to a
reference genome sequence in order for it to be usable by BWA. If you have created your own BWA index outside of DNAnexus,
place all the index files in a gzipped tar archive and provide that as the input; if you have a reference genome sequence in FASTA
format, you can index it with the BWA FASTA Indexer app). Some pre-indexed genomes are also available as suggested inputs. See
also ['which human reference sequence should I use?'](https://answers.dnanexus.com/p/183/) on DNAnexus Answers.

## What does this app output?

This app outputs the mappings, as a coordinate-sorted BAM file (`*.bam`). An associated BAM index file (`*.bai`) is also generated.

## How does this app work?

This app uses the BWA-MEM algorithm from the BWA software package. For general information, consult the BWA manual at:

http://bio-bwa.sourceforge.net/bwa.shtml

It also uses scripts developed by [Arima Genomics](http://www.arimagenomics.com/) to identify chimeric reads and merge
the bam files from mapping the forward and reverse reads.

This app performs the following steps:

- Mapping each of the forward and reverse reads separately with `bwa mem`.
- Conversion to BAM with `samtools view`.
- Filtering of each bam file for chimeric reads
- Merging of the two filtered bam files.
- Sorting, de-duplication, and index file generation using bamsormadup.
