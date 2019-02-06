# fArcCen1 Assembly Tutorial

This tutorial covers end-to-end assembly of the species [Flier Cyclid](https://vgp.github.io/genomeark/Archocentrus_centrarchus/]) `fArcCen1`. It covers 
some of the workflows and best practices for running your own VGP assembly.

## Getting Started

When assigned an assembly, you will given a species ID for your given genome 
(in this example `fArcCen1`) and a project will be shared with you containing 
the raw files for the genome which have been linked from the VGP AWS bucket.

The project will have the Species ID name (`fArcCen1`) and contain the following files:

```
.
└── genomic_data
   ├── phase (or arima)
   │   ├── fArcCen1_DDHiC_R1.fastq.gz
   │   ├── fArcCen1_DDHiC_R2.fastq.gz
   │   ├── fArcCen1_S3HiC_R1.fastq.gz
   │   └── fArcCen1_S3HiC_R2.fastq.gz
   ├── pacbio
   │   ├── m54062_170721_074508.subreads.bam
   │   ├── m54062_170721_175442.subreads.bam
   │   ├── m54062_170722_040411.subreads.bam
   │   ├── m54062_170809_134833.subreads.bam
   │   ├── m54062_170809_234737.subreads.bam
   │   ├── m54142_170706_115228.subreads.bam
   │   ├── m54142_170707_122532.subreads.bam
   │   ├── m54142_170711_151050.subreads.bam
   │   ├── m54142_170712_212934.subreads.bam
   │   ├── m54142_170810_135956.subreads.bam
   │   ├── m54142_170810_235924.subreads.bam
   │   ├── m54144_170808_141859.subreads.bam
   │   ├── m54144_170809_001755.subreads.bam
   │   ├── m54144_170811_141157.subreads.bam
   │   ├── m54144_170812_001052.subreads.bam
   │   ├── m54144_170817_111904.subreads.bam
   │   ├── m54144_170817_211802.subreads.bam
   │   └── m54144_170818_072740.subreads.bam
   ├── 10x
   │   ├── fArcCen1_S1_L001_I1_001.fastq.gz
   │   ├── fArcCen1_S1_L001_R1_001.fastq.gz
   │   ├── fArcCen1_S1_L001_R2_001.fastq.gz
   │   ├── fArcCen1_S1_L002_I1_001.fastq.gz
   │   ├── fArcCen1_S1_L002_R1_001.fastq.gz
   │   ├── fArcCen1_S1_L002_R2_001.fastq.gz
   │   ├── fArcCen1_S1_L003_I1_001.fastq.gz
   │   ├── fArcCen1_S1_L003_R1_001.fastq.gz
   │   ├── fArcCen1_S1_L003_R2_001.fastq.gz
   │   ├── fArcCen1_S1_L004_I1_001.fastq.gz
   │   ├── fArcCen1_S1_L004_R1_001.fastq.gz
   │   └── fArcCen1_S1_L004_R2_001.fastq.gz
   └── bionano
       ├── fArcCen1_Saphyr_DLE-1.cmap.gz
       ├── fArcCen1_Saphyr_DLE-1_1265239.bnx.gz
       └── fArcCen1_Saphyr_DLE-1_1265240.bnx.gz
```

This includes raw data from 4 data types:
1. HiC (provided by Phase Genomics or Arima) (`*.fastq.gz`)
2. Pacbio Sequel reads (`*.bam`)
3. 10X Genomics linked reads (`*.fastq.gz`)
4. Bionano optical maps (`*.cmap`)

To make sure the project is configured correctly, navigate to the `Settings` tab of the project. In addition to the project name (`fArcCen1`), the following fields should be configured as such:

```
Billed To: org-erich_lab
Region: AWS (US East)
Tags: White-faced whistling Duck, vgp, vgl
```

The `Tags` field will contain the common species name which is specific to your species.

Note: All work should be done in the project shared with you. Do **not** create a new project to work in.

## Falcon and Unzip Assembly

Navigate to the **VGP Tools** folder and locate the **Falcon and Unzip Workflow**. Right click the workflow and select `Copy`, after which a dialogue box will pop up requesting that you select a project. Select your working project (`fArcCen1`).

In your working project, click the workflow to open it in `Run` mode. Look through the workflow to make sure all instances and inputs are configured correctly.

The following are good to check as they tend to be misconfigured:
1. Under the "Unzip Track Reads" stage, the instance type should be set to `mem4_ssd1_x128`

Before configuring the workflow, it is good practice to create an edit-able copy of the workflow in case anything is misconfigured or needs to be re-run. This can be done by selecting `Workflow Actions` from the `Run Analysis...` screen and selecting the `Save Template` action. This will take you to a copy of the workflow that can be modified. Workflows are automatically saved if any changes are made. Return the the `Run Analysis...` screen by clicking the `Start Analysis` button.

![](http://)

Once the workflow is configured, select the `BAM Files` input under the `BAM to FASTA` stage. This will pop up a dialogue to select input files. Select the *PacBio Sequel Files* as input.

Under the `Pre-Assembly Calculate Read Length Distribution` stage, click the gear icon to open the parameters panel and fill in the `Genome Size` parameter with the given species' expected genome size. For `fArcCen1`, the estimated genome size is 0.99Gbp. This can also be estimated by running `Jellyfish and Genomescope` on the 10x Genomics reads.

In addition to selecting parameters, it is useful to specify an output folder for the workflow. Under `Workflow Actions`, select `Set Output Folder`. Create a new folder with the name `assembly_version1.5` and select this as the output folder for the `Falcon and Unzip Assembly Workflow`. All of the output folders for the individual stages should already be configured.

The workflow should now be in the "Runnable" state. Before running, make sure to save your workflow changes (including input specification) by selecting `Workflow Actions` and selecting `Update Workflow With Changes`. Finally, click `Run as Analysis...` to launch the workflow.

Should any failures occur, the workflow can be restarted from the point of failure by navigating to the "Monitor" tab of the project and selecting the launched Falcon and Unzip Assembly analysis. There you can select the "Rerun As Analysis" button. Be sure to wait while all inputs and links populate before making any changes to the analysis parameters.

## Scaffolding

### 1. Purge Haplotigs

The Purge Haplotigs step consists of two stages:
1. Stage 1 performs alignment with Minimap2 and generates a coverage histogram plot
2. Stage 2 uses cutoffs to purge haplotigs

From the `VGP Tools` project, copy the `Purge Haplotigs Stage 1` and `Purge Haplotigs Stage 2` workflows. Select the `Stage 1` workflow.

Under the `Rename Contigs` stage, select the `*.fasta` input file and locate the assembly `*.fasta` file located under `Unzip Stage 5`.

Under the `Minimap2` stage, select the *PacBio Sequel Files* as input.