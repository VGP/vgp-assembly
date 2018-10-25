<!-- dx-header -->
# Salsa (DNAnexus Platform App)

SALSA: A tool to scaffold long read assemblies with Hi-C data (https://github.com/machinegun/SALSA)

## What Does this App Do?

This app will run checkout/compile SALSA and scaffold provided contigs using HiC data.

## How Does this App Work?

This app uses the SALSA with an input contig assembly and input bam files. It breaks errors in contigs and scaffolds them based on the links found in the bam file. An optional GFA can be provided if it is available for the assembly. The names of the sequences in the GFA file must match the names in the fasta file. Output from all iterations (fasta and agp) is saved.
