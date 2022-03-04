[decontam_pipe.pdf](https://github.com/Nadolina/decontam_pipe_2/files/8181224/decontam_pipe.pdf)

Following the VGP genome assembly pipeline 2.1, scaffolded assemblies require decontamination to remove any terminal gaps, non-target contaminants (i.e./ bacterial, human contamination) and sequences originating from the mitochondria. This pipeline has been tested on vertebrate genomes only. 

**Dependencies**

The pipe is written in shell and python scripts.
- Python3
- Python modules Biopython and Pandas 
- dustmasker v1.0.0
- Kraken2; standard plus database (https://github.com/DerrickWood/kraken2/wiki/Manual#kraken-2-databases)
- blast+/blastn v2.5.0+
- NCBI mitochondria database (https://ftp.ncbi.nlm.nih.gov/refseq/release/mitochondrion/mitochondrion.2.1.genomic.fna.gz)

**Inputs and running the pipeline**

All scripts/programs are submitted via the shell pipeline. Note that the pipe was developed for internak VGL/VGP use and so is based on submission to a slurm queueing system. The pipe can be run from the commandline:

sbatch -p <parition> VGP_decontamination_pipe.sh <fasta> <unique ID>
  
Fasta files must be decompressed. The unique ID can be anything, it is just used for naming throughout the pipe; for VGP purposes, we use the TOLID (i.e./ bTaeGut2).

**Outputs**
1. class-/ unclass_bTaeGut2_<fastaname> - two separate files containing the classified (contaminant) and unclassified (target) scaffolds; a future update will include the removal of the unclassified file since a final fasta is generated at the end of the pipe.
2. contam_scaffs_<unique ID>.txt - compiled list of contaminant scaffolds from the kraken2 and mito-blast subprocesses 
3. mito_blast_<fastaname>.report - from the parse_mito_blast.py script which summarizes the blast output table, listing the highest coverage scaffold-accession number pairs (high coverage = mito-contaminant)
4. N_sub_masked_<fastaname> - fasta with hard masking after dustmasker + sub_soft_hard_mask.py 
5. trimmed_<fastaname> - THE FINAL OUTPUT; a scaffolded assembly from which terminal gaps and all contaminant (non-target and mito) have been removed 

