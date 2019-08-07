# Minimap2 Long Read Aligner

## What does this app do?

This app maps PacBio or Oxford Nanopore long reads to a reference genome
with the Minimap2 aligner.

## What are typical use cases for this app?

This app can be used to map reads to a reference genome, which is a typical step in most bioinformatics
analyses. It is suitable as a first step if you are planning on performing assembly polish or evaluation via mappings.

This app uses Minimap2, a versatile  tool for aligning sequencing reads to long reference sequences.

## What data are required for this app to run?

This app requires reads files in one of the following formats:
* PacBio files in `*.subreads.bam` format
* PacBio files in `*.fasta.gz` or `*.fastq.gz` format
* ONT files in `*.fasta.gz` or `*.fastq.gz` format

The app also requires a reference genome sequence in FASTA format or as a pre-compiled
Minimap2 index in `*.mmi` format. If a `*.mmi` file is not provided, one is generated and output by the app.

In addition, users can specify which pre-set mapping parameters to use by 
selecting the reads datatype. The `map-pb` preset is used for PacBio reads and
the `map-ont` preset is used for Oxford Nanopore reads.

## What does this app output?

This app outputs the mappings, as a coordinate-sorted BAM file (`*.bam`), and 
a corresponding BAM mappings index (`*.bai`) file.

## How does this app work?

This app uses Minimap2. For general information, consult the Minimap2 manual at:

https://lh3.github.io/minimap2/

This app performs the following steps:

- Mapping with `minimap2`
- Sorting, indexing and marking duplicates with `biobambam2`
- If the `Run pbbamify` is enable, the pbbamify will be run after pbmm2 to enable compatibility with Arrow. The current pbmm2 from Conda does not provide md5 for contigs which is required for polishing using standard PacBio release code.
