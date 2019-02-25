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

Before configuring the workflow, it is good practice to create an edit-able copy of the workflow in case anything is 
misconfigured or needs to be re-run. This can be done by selecting `Workflow Actions` from the `Run Analysis...` screen 
and selecting the `Save Template` action. This will take you to a copy of the workflow that can be modified. Workflows 
are automatically saved if any changes are made. Return the the `Run Analysis...` screen by clicking the `Start Analysis` 
button.

![Workflow Actions: Save as Template](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/WorkflowSaveAsTemplate.png)

Once the workflow is configured, select the `BAM Files` input under the `BAM to FASTA` stage. This will pop up a dialogue to select input files. Select the *PacBio Sequel Files* as input.

![BAM to FASTA input](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/BAMtoFASTAinput.png)

Under the `Create Raw Reads Dazzler DB` stage, click the gear icon to open the parameters panel and fill in the 
`Genome Size` parameter with the given species' expected genome size. For `fArcCen1`, the estimated genome size is 
0.99Gbp, so we fill in `0.99G`. This can also be estimated by running `Jellyfish and Genomescope` on the raw 10x 
Genomics reads.

![Create Raw Reads stage](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/CreateRawReadsConfig.png)
![Configure genome size](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/LengthCutoffConfig.png)

In addition to selecting parameters, it is useful to specify an output folder for the workflow. Under `Workflow Actions`, 
select `Set Output Folder`. Create a new folder with the name `assembly_version1.5` and select this as the output folder 
for the `Falcon and Unzip Assembly Workflow`. All of the output folders for the individual stages should already be 
configured. When the analysis is complete, there will be a total of 10 subdirecories in your specified output directory, 
one for each stage of the workflow:

1. BAM to FASTA (bam_to_fasta)
2. FALCON stage 0 (stage_0)
3. FALCON stage 1 (stage_1)
4. FALCON stage 2 (stage_2)
5. FALCON stage 3 (stage_3)
6. FALCON Unzip stage 1 (unzip_stage_1)
7. FALCON Unzip stage 2 (unzip_stage_2)
8. FALCON Unzip stage 3 (unzip_stage_3)
9. FALCON Unzip stage 4 (unzip_stage_4)
10. FALCON Unzip stage 5 (unzip_stage_5)

The workflow should now be in the "Runnable" state. Before running, make sure to save your workflow changes (including 
input specification) by selecting `Workflow Actions` and selecting `Update Workflow With Changes`. This will make it 
 easier to modify and relaunch the workflow should any failures occur. Finally, click `Run as Analysis...` to launch 
 the workflow.

![Remember to save your workflow configuration by clicking "Update Workflow with Changes"](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/WorkflowUpdateWithChanges.png)
![Click the Run as Analysis button to launch the workflow.](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/RunAsAnalysisGo.png)

Should any failures occur, the workflow can be restarted from the point of failure by navigating to the "Monitor" tab of
 the project and selecting the launched Falcon and Unzip Assembly analysis. There you can select the "Rerun As Analysis" 
 button. Be sure to wait while all inputs and links populate before making any changes to the analysis parameters.

The final output should look like this:
```
assembly_v1.5
├── bam_to_fasta
│   ├── m#####_######_######.subreads.fasta.gz
│   ├── m#####_######_######.subreads.fasta.gz
│   ├── m#####_######_######.subreads.fasta.gz
│   └── ...
├── stage_0
│   ├── read_counts.csv
│   └── read_length_distribution.pdf
├── stage_1
│   ├── genome_scope
│   │   ├── mer_counts.csv
│   │   ├── mer_counts.gs.log.png
│   │   ├── mer_counts.gs.png
│   │   ├── mer_counts.model.txt
│   │   ├── mer_counts.progress.txt
│   │   └── mer_counts.summary.txt
│   ├── out.####.fasta.gz
│   ├── out.####.fasta.gz
│   ├── out.####.fasta.gz
│   ├── out.####.fasta.gz
│   ├── out.####.fasta.gz
│   ├── ...
│   ├── preads_rep.tar.gz
│   ├── preads_tan.tar.gz
│   ├── raw_reads.###.las
│   ├── raw_reads.###.las
│   ├── raw_reads.###.las
│   ├── raw_reads.###.las
│   ├── ...
│   ├── raw_reads_rep.tar.gz
│   ├── raw_reads_tan.tar.gz
│   ├── read_counts.csv
│   └── read_length_distribution.pdf
├── stage_2
│   ├── preads.###.las
│   ├── preads.###.las
│   ├── preads.###.las
│   └── ...
├── stage_3
│   ├── a_ctg.fasta.gz
│   ├── a_ctg_tiling_path.gz
│   ├── asm.gfa.gz
│   ├── asm.gfa2.gz
│   ├── contig.gfa2.gz
│   ├── ctg_paths.gz
│   ├── p_ctg.fasta.gz
│   ├── p_ctg_tiling_path.gz
│   ├── preads4falcon.fasta.gz
│   ├── read_counts.csv
│   ├── read_length_distribution.pdf
│   ├── sg.gfa.gz
│   ├── sg.gfa2.gz
│   ├── sg_edges_list.gz
│   └── utg_data.gz
├── unzip_stage_1
│   ├── ctg_list
│   ├── pread_ids.gz
│   ├── pread_to_contigs.gz
│   ├── rawread_ids.gz
│   ├── rawread_to_contigs.gz
│   ├── read_to_contig_map.gz
│   ├── reads_fastas_##.tar.gz
│   ├── reads_fastas_##.tar.gz
│   ├── reads_fastas_##.tar.gz
│   ├── reads_fastas_##.tar.gz
│   ├── ...
│   ├── ref_fastas_###.tar.gz
│   ├── ref_fastas_###.tar.gz
│   ├── ref_fastas_###.tar.gz
│   └── ...
├── unzip_stage_2
│   ├── all_phased_reads.gz
│   ├── rid_to_phase.###.tar.gz
│   ├── rid_to_phase.###.tar.gz
│   ├── rid_to_phase.###.tar.gz
│   ├── rid_to_phase.###.tar.gz
│   ├── rid_to_phase.###.tar.gz
│   ├── rid_to_phase.###.tar.gz
|   └── ...
├── unzip_stage_3
│   ├── all_h_ctg.fa.gz
│   ├── all_h_ctg_edges.gz
│   ├── all_h_ctg_ids.gz
│   ├── all_p_ctg.fa.gz
│   ├── all_p_ctg_edges.gz
│   ├── asm.gfa.gz
│   └── sg.gfa.gz
├── unzip_stage_4
│   ├── contig_groups.json
│   ├── mappedread_###.tar.gz
│   ├── mappedread_###.tar.gz
│   ├── mappedread_###.tar.gz
│   ├── mappedread_###.tar.gz
│   └── ...
└── unzip_stage_5
    ├── cns_h_ctg.fasta.gz
    ├── cns_h_ctg.fastq.gz
    ├── cns_p_ctg.fasta.gz
    └── cns_p_ctg.fastq.gz

```

The "unzip_stage_5" directory contains the primary assembly (`cns_p_ctg.fasta.gz`) and the haplotype assembly 
(`cns_h_ctg.fasta.gz`) which will be used in the scaffolding workflow.

## Scaffolding

### 1. Purge Haplotigs

The Purge Haplotigs step consists of two stages:
1. Part 1 performs alignment with Minimap2 and generates a coverage histogram plot
2. Part 2 uses cutoffs to purge haplotigs

From the `VGP Tools` project, copy the `Purge Haplotigs Part 1` and `Purge Haplotigs Part 2` workflows to your working 
project. In the `fArcCen1` working project, select the `Purge Haplotigs Part 1` workflow.

Under the `Rename Contigs` stage, select the `*.fasta` input file by navigating to `assembly_v1.5/unzip_stage_5` and 
selecting the file `cns_p_ctg.fasta.gz`. This represents the primary unzipped assembly.

Under the `Minimap2` stage, provide the *raw* PacBio FASTA files as input by navigating to `assembly_v1.5/bam_to_fasta` 
and selecting the full list of FASTA files.

Finally, under "Workflow Actions", select "Set Output Folder" and create a new folder for "Scaffolding", selecting the 
`Scaffolding` folder as the output folder. All of the apps should have individual output folders configured. Now run the
`Purge Haplotigs Stage 1` workflow.

Once the analysis completes, take a look at the `cns_p_ctg.renamed.coverage_histo.pdf` output file by selecting the file
and clicking on "Preview".

![Right click the output file and select "Preview"]()

A preview window will show you a PDF containing a bimodal histogram which can be used to select the Purge Haplotigs cutoffs.
If the histogram is not biomodal (contains only one peak), you can skip ahead to the Scaff10x step.

If the histogram is bimodal, select the cutoffs as follows:
* _Low_: Threshold before the first peak
* _Mid_: Threshold between two two peaks
* _High_: Threshold above the second peak

![Coverage histogram with indicated cutoffs](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/PurgeHaplotigsReadhist.png)

Once you've determined the appropriate cutoffs, select the `Purge Haplotigs Part 2` workflow. Specify the file inputs as
follows:
* `Reference Genome`: "renamed" reference genome from the "scaffolding/purge_haplotigs" directory (`cns_p_ctg.renamed.fasta.gz`)
* `Alignment`: merged mapping output file from the "scaffolding/purge_haplotigs" directory (`m54.bam`)
* `Genome Coverage`: The `*.gencov` output file from the "scaffolding/purge_haplotigs" directory (`cns_p_ctg.renamed.gencov`)

In addition, click the gear icon to the right of the app to open the parameters panel and specify the cutoffs. For 
`fArcCen1`, we've chosen the following cutoffs: Low = 10, Mid = 48, High = 185.

![Cutoffs specified in Part 2 app](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/PurgeHaplotigsCutoffs.png)

Once the inputs and parameters are specified, specify the Output Folder as before to the "Scaffolding" folder. Now run 
the analysis.

The final output should look like this:
```
scaffolding
└── purge_haplotigs
    ├── mapped_reads
    │   ├── cns_p_ctg.renamed.mmi
    │   ├── m#####_######_######.subreads.mapped.bam
    │   ├── m#####_######_######.subreads.mapped.bam.bai
    │   ├── ...
    │   └── ...
    ├── cns_p_ctg.renamed.coverage_histo.pdf
    ├── cns_p_ctg.renamed.fasta.gz
    ├── cns_p_ctg.renamed.gencov
    ├── cns_p_ctg.renamed_coverage_stats.csv.gz
    ├── curated.artefacts.fasta.gz
    ├── curated.contig_associations.log
    ├── curated.fasta.gz
    ├── curated.haplotigs.fasta.gz
    ├── curated.reassignments.tsv
    ├── dotplots_reassigned_contigs.tar.gz
    ├── dotplots_unassigned_contigs.tar.gz
    ├── m54.bam
    └── m54.bam.bai
```

In this output, the primary contigs are contained in `curated.fasta.gz` and the purged haplotigs are contained in 
`curated.haplotigs.fasta.gz`.

## 2. Scaff10x Workflow

The Scaff10x workflow performs two rounds of scaffolding with Scaff10x using the 10X Genomics raw reads. Copy the 
`Scaff10x 2 Rounds` workflow from the `VGP Tools` directory.

Select the file inputs for the first step of the workflow:
* assemble_genome_fastagz: This should be the _purged_ primary contigs from Purge Haplotigs (`curated.fasta.gz`)
* R1 reads: Under `genomic_data/10x`, all reads containing `_R1_*.fastq`
* R2 reads: Under `genomic_data/10x`, all reads containing `_R2_*.fastq`

Under the first round parameters, make sure the `is_raw` parameter is specified to `True`. Under both round 1 and round 
2 parameters you may also specify the output prefix for the output files if desired.

In addition, specify the `Output Folder` for the workflow to `scaffolding` as before. Launch the analysis.

![Configured Scaff10x workflow](https://raw.githubusercontent.com/VGP/vgp-assembly/master/tutorials/images/Scaff10XRunnable.png)

The final output of scaff10x should look like this:
```
scaffolding/scaff10x/
├── round1
│   ├── read-BC_1.fastq.gz
│   ├── read-BC_2.fastq.gz
│   └── scaffolds.fasta.gz
└── round2
    ├── read-BC_1.fastq.gz
    ├── read-BC_2.fastq.gz
    └── scaffolds.fasta.gz
```

The output under `round2/scaffolds.fasta.gz` corresponds to `fArcCen1_s1.fasta.gz` using the VGP naming convention.

## 3. Bionano Hybrid Scaffolding

In this step we will be using the Bionano assembled CMAP data together with the scaffold `fArcCen1_s1.fasta.gz` from 
step 2 to perform hybrid scaffolding on the primary haplotig.

Before copying the workflow, we need to first take a look at the Bionano input data:
```
genomic_data/bionano/
├── fArcCen1_Saphyr_DLE-1.cmap.gz
├── fArcCen1_Saphyr_DLE-1_1265239.bnx.gz
└── fArcCen1_Saphyr_DLE-1_1265240.bnx.gz
```

From the input data, we can see that there is a single `*.cmap.gz` file generated using the `DLE-1` enzyme. Therefore,
copy the `Step 3 Bionano Hybrid Scaffolding 1 Enzyme` workflow from VGP tools. In some cases you may see two `*.cmap.gz`
files in your input data corresponding to two enzymes used for the Bionano data, and therefore will need to use the 2 enzyme workflow.

For the inputs select as follows:
* CMAP input: `fArcCen1_Saphyr_DLE-1.cmap.gz` under `genomic_data/bionano/`
* FASTA input: `scaffolds.fasta.gz` under `scaffolding/scaff10x/round2`

As before, specify the workflow output folder to `scaffolding`. A `bionano` folder will be automatically created under that.

The Bionano tool is prone to hanging if the memory requirements are not met. Therefore, make sure it is running on a 
large memory instance, such as the `mem3_ssd1_x32` instance. The app is configured to aggressively time out if it does 
not complete in under 6 hours. If this happens, rerun the analysis on a larger memory instance.

Once the app completes, the output should look as follows:

As part of the Bionano hybrid scaffolding, you should also perform the following post-processing steps:
1. Concatenate `UNSCAFFOLDED` and `SCAFFOLD_FINAL_NCBI` outputs
2. Remove leading and trailing N's from scaffolded output.

To concatenate the two outputs, select the "Run Analysis" button in your project and select the "File Concatenator" app.

Supply the inputs `UNSCAFFOLDED` and `SCAFFOLD_FINAL_NCBI` to the app as well as the output name `fArcCen1_s2.fasta.gz`.

## Step 4.