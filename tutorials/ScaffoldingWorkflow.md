# VGP Scaffolding Workflow
## Version: 1.5

This workflow is run after Falcon and Falcon Unzip assembly to do scaffolding with Hi-C, Bionano and 10X Genomics data types.

The outputs of Falcon and Falcon Unzip are:
* Primary contigs assembly (p1.fasta.gz)
* Haplotype contigs assembly (q2.fasta.gz)

### 1. Purge Haplotigs

Purge Haplotigs [https://bitbucket.org/mroachawri/purge_haplotigs] identifies heterozygous regions which were mistaken for primary contigs instead of alternate and removes them.

*INPUT*: Primary contig assembly (p1.fasta.gz), raw PacBio reads.
*Steps*:
0. Remove '|' characters from Unzip contigs (these are added by Arrow and break Purge Haplotigs)
1. Align PacBio reads to Primary Contig using *Minimap2*
2. Merge output BAM files into a single Merged BAM
3. Create histogram of read coverage using *purge_haplotigs readhist*
3. Manually determine threshold by looking at histogram and finding where the second peak cutoffs are (If there's no second peak, you can continue to next scaffolding step.)
4. Purge haplotigs by running *contigcov_and_purge*

*OUTPUT*: Purged p1.fasta.gz

DNAnexus Apps:
* Rename Contigs [https://platform.dnanexus.com/app/vgp_rename_contigs]
* Minimap2 Long Read aligner [https://platform.dnanexus.com/app/minimap2_align_longread]
* Bamtools Mappings Merge [https://platform.dnanexus.com/app/bamtools_merge]
* Purge Haplotigs Readhist [https://platform.dnanexus.com/app/purge_haplotigs_readhist]
* Purge Haplotigs Contigcov and Purge [https://platform.dnanexus.com/app/purge_haplotigs_contigcov_and_purge]

### 2. scaff10x:

Scaff10x [https://github.com/wtsi-hpag/Scaff10X] uses 10X genomic data to scaffold the primary assembly.
Use the scaff10x workflow under [VGP Tools] project to do 2-rounds of scaffolding.

*Input Data*: Purge primary contig assembly, 10x Genomics raw reads
1. Run 1 round of scaff10x using `"raw_input": true` parameter
2. Use outputs from 1st round to run 2nd round of scaff10x using `"raw_input": false` parameter

*OUTPUT*: Scaffolded FASTA (s1.fasta.gz)

*Note:* Break10X is not used.

DNAnexus Apps:
* Scaff10x [https://platform.dnanexus.com/app/scaff10x]

### 3. Bionano hybrid scaffolding:

Bionano hybrid scaffolding uses Bionano optical maps to scaffold the primary assembly. An *assembled optical map* (`*.cmap`) is required to perform hybrid scaffolding.

*INPUT*: Scaffold from Scaff10x (s1.fasta.gz), Assembled Bionano CMAP (BspQI + BssSI or DLE1)
You will also need to know which enzymes the optical maps were generated with.
To run Bionano hybrid scaffolding, select the 1 or 2 enzyme app. Use 1 enzyme for DLE1 optical maps or for either BspQI or BssSI optical maps. For the 2 enzyme optical map protocol, use the 2 enzyme hybrid scaffolding app.

1. Run the appropriate hybrid scaffolding app
2. Concatenate [https://platform.dnanexus.com/app/file_concatenator] the NCBI_SCAFFOLDED.fasta output with the UNSCAFFOLDED.fasta output to create s2.fasta.gz
3. Trim off the leading / trailing N's in s2.fasta.gz

OUTPUT: s2.fasta.gz

DNAnexus Apps:
* Bionano 1 Enzyme Hybrid Scaffolding [https://platform.dnanexus.com/app/bionano_hybrid_1enzyme]
* Bionano 2 Enzyme Hybrid Scaffolding [https://platform.dnanexus.com/app/bionano_hybrid_2enzyme]

### 3. Salsa scaffolding:

Salsa scaffolding uses Hi-C data to scaffold the assembly.
Use the Arima mapping pipeline [https://github.com/ArimaGenomics/mapping_pipeline] and Salsa2.1 [https://github.com/machinegun/SALSA/releases/tag/v2.1]

INPUT: s2.fasta.gz and HiC raw reads (`*.fastq`). You will also need the restriction sequence (such as `GATC`) used to generate the HiC reads.

1. Align the HiC reads using the Arima mapping pipeline
2. Run Salsa2 on aligned reads and the s2.fasta.gz scaffold

OUTPUT: Salsa scaffold (s3.fasta.gz)

3. Concatenate the Salsa2 output with the Haplotigs (q2.fasta.gz) from FALCON Unzip

OUTPUT: Salsa scaffold + haplotigs (s4.fasta.gz)

DNAnexus Apps:
* Arima Mapping [https://platform.dnanexus.com/app/arima_mapping]
* Salsa 2 [https://platform.dnanexus.com/app/salsa]
* File Concatenator [https://platform.dnanexus.com/app/file_concatenator]

### 5. Arrow polishing:

Polish the scaffold using PacBio raw reads and the *Arrow* algorithm

* Use arrow with blasr mapping for genomes <4Gb (https://github.com/skoren/ArrowGrid/) 
* Use minimap2 arrow for genomes >4Gb (https://github.com/skoren/ArrowGrid/tree/no_blasr)

INPUT: Salsa scaffold with Haplotigs (s4.fasta.gz) and PacBio raw reads (`*.subreads.bam` or `*.subreads.fasta`)
1. Align the raw reads to the scaffold using PBalign (Blasr) or Minimap2
2. If aligned using Minimap2, run PBbamify (https://github.com/PacificBiosciences/pbbam/wiki/pbbamify) to make aligned BAM files compatible with Arrow.
3. Polish the reads using PacBio's Arrow algorithm

OUTPUT: Arrow consensus FASTA (t1.fasta.gz)

DNAnexus Apps:
* Minimap2 long read aligner [https://platform.dnanexus.com/app/minimap2_align_longread]
* BLASR aligner [https://platform.dnanexus.com/app/pbalign]
* Arrow Polish [https://platform.dnanexus.com/app/run_polish]

### 6. Freebayes polishing:

Use 10X Genomics raw reads to polish the assembly.
10X Longranger [https://support.10xgenomics.com/genome-exome/software/overview/welcome] is used for alignment, Freebayes variant caller [https://github.com/ekg/freebayes] is used for variant detection and BCFtools [https://samtools.github.io/bcftools/bcftools.html] is used for calling consensus.

*Note:* Run this workflow two times
*INPUT*: t1.fasta.gz and 10X Genomics Raw Reads (`*.fastq`)

1. Create reference genome for 10X Longranger aligner using `longranger mkref`
2. Align 10x Genomics `*.fastq` files to reference using `longranger align`
3. Call variants using `freebayes`
4. Create consensus FASTA using `bcftools consensus`
5. Run the genome coverage step in the workflow at the end, so we can calculate the QV.
6. (Second round) Repeat steps 1 - 5 for the second round using t2.fasta.gz as input.

OUTPUT1: t2.fasta.gz (1st round)
OUTPUT2: t3.fasta.gz (2nd round)
and QVs for each previous steps

DNAnexus Apps:
* Longranger mkref [https://platform.dnanexus.com/app/10x_longranger_mkref]
* Longranger align [https://platform.dnanexus.com/app/10x_longranger_align]
* Freebayes Variant Caller [https://platform.dnanexus.com/app/freebayes]
* BCFtools Consensus [https://platform.dnanexus.com/app/bcftools_consensus]

### Data Transfer
All OUTPUT and INTERMEDIATE files needs to be transferred back to AWS.

The OUTPUT goes under s3://genomeark/species/species_name/species_id/assembly_v1.5/ .
The INTERMEDIATES go under s3://genomeark/species/species_name/species_id/assembly_v1.5/intermediates/

DNAnexus Apps:
* DNAnexus to VGP S3 Exporter [https://platform.dnanexus.com/app/dx_to_vgp_s3_file_transfer]

### Future Work

Proposed streamlined workflow:
1. Single workflow runs Falcon + Falcon Unzip + Purge Haplotigs Step 1
2. Second Workflow runs Purge Haplotigs through Freebayes Polish

Some choices would be removed from the workflow:
* Bionano hybrid scaffolding will always be 1 enzyme (DLE)
* Minimap2 will be used for Arrow polishing for all asssemblies


### Notes

Apps are available to `org-vgp` members on DNAnexus. Please contact VGP administrators to be added.
