# _fArcCen1_ Assembly Tutorial

This tutorial covers the assembly of the fish species [Flier Cyclid](https://vgp.github.io/genomeark/Archocentrus_centrarchus) (_Archocentrus centrarchus_) using the [DNAnexus platform](https://platform.dnanexus.com/).
The overall assembly pipeline can be depicted in the following simplified diagram:

![DNAnexus assembly diagram 1.6](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/DNAnexus_VGP_1.6_diagram_28012020.png)


**Index:**

. [Getting Started](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#getting-started)

. [Falcon and Unzip Assembly](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#1-falcon-and-unzip)

. [Purge Dups](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#2-purge-dups)

. [10X Scaffolding](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#1-10x)

. [Bionano Hybrid Scaffolding](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#2-bionano)

. [Salsa Scaffolding](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#3-salsa)

. [Arrow Polishing](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#1-arrow)

. [Freebayes Polishing](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#2-freebayes)

. [Final Checks](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#final-checks)


IMPORTANT: Remember that your questions are always welcome in the ["training" channel of Slack](https://genomeark.slack.com/archives/CE7FU8YAC)!

<br/>

## Getting Started

When assigned an assembly, you will be given a **Species ID** for your given genome (in this example _**fArcCen1**_) and a project will be shared with you containing the raw files for the genome which have been linked from the VGP AWS bucket.

The root folder of the project will have the _Species ID_ name (_fArcCen1_) and contain the following folders and files:

```
fArcCen1
└── genomic_data
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
    ├── bionano
    │   ├── fArcCen1_Saphyr_DLE-1.cmap.gz
    │   ├── fArcCen1_Saphyr_DLE-1_1265239.bnx.gz
    │   └── fArcCen1_Saphyr_DLE-1_1265240.bnx.gz
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
    └── phase (or arima)
        ├── fArcCen1_DDHiC_R1.fastq.gz
        ├── fArcCen1_DDHiC_R2.fastq.gz
        ├── fArcCen1_S3HiC_R1.fastq.gz
        ├── fArcCen1_S3HiC_R2.fastq.gz        
        └── re_bases.txt
```

This includes 4 types of **raw** data:
1. 10X Genomics linked reads (`*.fastq.gz`). The reads are contained in the compressed _fastq_ files _\_R1\__ and _\_R2\__ (the files _\_I1\__ are indexes, not to be used). Shortly, 10X reads are Illumina short reads that contain a barcode that link each one of them to the DNA molecule (ideally a chromosome) from where they come (for more details, click [here](https://www.10xgenomics.com/linked-reads/)). 
2. Bionano optical maps (`*.cmap`). 
3. Pacbio Sequel reads (`*.bam`). Binary files that contain long reads generated with the Sequel platform (for more information, click [here](https://www.pacb.com/products-and-services/sequel-system/)). Please be sure to only use the _.subreads_ files.
4. HiC (provided by Arima or Phase Genomics) (`*.fastq.gz`).

<br/>

To make sure the project is configured correctly, navigate to the "Settings" tab of the project. In addition to the project name (_**fArcCen1**_), the following fields should be configured as such:

```
Billed To: org-erich_lab
Region: AWS (US East)
Tags: Flier Cyclid
```


IMPORTANT: All work should be done in the project shared with you. Do **not** create a new project to work in.

<br/>

## Assembly

### 1. Falcon and Unzip

FALCON and FALCON-Unzip are _de novo_ genome assemblers for PacBio long reads (for more information, click [here](https://pb-falcon.readthedocs.io/en/latest/about.html#overview)).

Since the Falcon workflow needs an estimated genome size, and it is useful to know some other properties of the genome before starting, you must first run the **meryl+genomescope_10x** workflow.

In your working project, click the green button `+ Add Data` and search and select **VGP Tools** in the "Other Project" tab. Search and select the latest version of the **meryl+genomescope_10x** workflow and click the green button `Add Data`, after which a dialogue box will pop up with a progress bar indicating that the workflow has been copied to the current location of your working project (the latest version of the workflows and applets should allways be in the main _VGP tools_ folder, make sure not use archived versions).

![DNAnexus working project](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/DNA_nexus_front.png)

<br/>

Click the workflow to open it in _Run_ mode. The workflow only needs to be configured for input files in the `Remove gembarcodes from 10x reads` stage. To add the input files, select all the `fastq.gz` in the `10x` folder (only _R1_ or _R2_ files). Next, to specify an output folder for the workflow, under `Workflow Actions`, select `Set Output Folder`, and create a folder named `assembly_vgp_standard_1.6`. Inside that new folder, create a folder named `meryl_genomescope`. Finally, click `Run as Analysis...` to launch the workflow.

At this stage, the final folder structure should look in general like this:
```
fArcCen1
├── assembly_vgp_standard_1.6
│   └── meryl_genomescope
│       └── ...
└── genomic_data
    └── ...
```

**!)** In addition to the genome size, it is useful to take a look to other values in the GenomeScope results plot (e.g. heterozygosity percentage). To further details about the GenomeScope plots and its interpretation, click [here](https://github.com/VGP/vgp-assembly/blob/master/tutorials/docs_1.6/Genomescope_overview_1.6.md#genomescope-plots-overview).

<br/>

Now, to formally start with the assembly workflow, click the green button `+ Add Data` and search and select **VGP Tools** in the "Other Project" tab. Search and select the latest version of the **vgp_falcon_and_unzip_assembly_workflow** and click the green button `Add Data`, after which a dialogue box will pop up with a progress bar indicating that the workflow has been copied to the current location of your working project.

In your working project, click the workflow to open it in _Run_ mode.

![Falcon and unzip workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/falcon_unzip_workflow.png)

<br/>

Before configuring the workflow, it is good practice to create an editable copy of the workflow in case anything is misconfigured or needs to be re-run. This can be done by selecting `Workflow Actions` and then selecting the `Save as template` action. This will take you to a copy of the workflow that can be modified. Workflows are automatically saved if any changes are made. Return to the `Run Analysis...` screen by clicking the `Start Analysis` button.

![Workflow Actions: Save as Template](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/WorkflowSaveAsTemplate.png)

<br/>

Look through the workflow to make sure all instances and inputs are configured correctly. Please check the following as it tend to be misconfigured: Under the `Unzip Track Reads` stage, the instance type should be set to `mem4_ssd1_x128`, unless something different is told to you in the training channel of Slack (for more information about memory and other instances in DNAnexus, click [here](https://github.com/VGP/vgp-assembly/blob/master/tutorials/docs_1.6/DNAnexus_instances.md)).

Once the workflow is configured, select the `BAM Files` input under the `BAM to FASTA` stage. This will pop up a dialogue window to select input files. Select the _PacBio Sequel Reads_ from the `pacbio` folder as input (as a good practice, please always select the corresponding files by locating them in their respective folders).

![BAM to FASTA input](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/BAMtoFASTAinput.png)

<br/>

Under the `Create Raw Reads Dazzler DB` stage, click the gear icon to open the parameters panel and fill in the "Estimated genome size" parameter with the given species' expected genome size. For `fArcCen1`, the estimated genome size is 0.99Gbp, so we fill in `0.99G`. The genome size can be obtained from the [Animal Genome Size Database](http://www.genomesize.com/) or if the species is no available there, can also be estimated using the result of the **meryl+genomescope_10x** workflow on the 10x Genomics reads as explained before).

![Create Raw Reads stage](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/CreateRawReadsConfig.png)
![Configure genome size](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/LengthCutoffConfig.png)

<br/>

In addition to selecting parameters, you should specify an output folder for the workflow. Under `Workflow Actions`, select `Set Output Folder`. Navigate to the folder `assembly_vgp_standard_1.6` and create a new folder named `intermediates`. Inside that new folder, create a folder named `falcon_unzip`. Select `falcon_unzip` as the output folder for the **vgp_falcon_and_unzip_assembly_workflow**. All of the output folders for the individual stages should already be configured. When the analysis is complete, there will be a total of 10 subfolders in your specified output folder, one for each stage of the workflow:

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

All the stages of the workflow should now be in the "Runnable" state. Before running, make sure to save your workflow changes (including input specification) by selecting `Workflow Actions` and selecting `Update workflow with changes`. This will make it easier to modify and relaunch the workflow should any failures occur. Finally, click `Run as Analysis...` to launch the workflow.

![Remember to save your workflow configuration by clicking "Update Workflow with Changes"](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/WorkflowUpdateWithChanges.png)
![Click the Run as Analysis button to launch the workflow.](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/RunAsAnalysisGo.png)

Should any failures occur, the workflow can be restarted from the point of failure by navigating to the "Monitor" tab of the project and selecting the launched **vgp_falcon_and_unzip_assembly_workflow**. There you can select the `Rerun As Analysis` button. Be sure to wait while all inputs and links populate before making any changes to the analysis parameters.

The final output should look like this:
```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       └── falcon_unzip
│           ├── bam_to_fasta
│           │   ├── m#####_######_######.subreads.fasta.gz
│           │   ├── m#####_######_######.subreads.fasta.gz
│           │   └── ...
│           ├── stage_0
│           │   ├── read_counts.csv
│           │   └── read_length_distribution.pdf
│           ├── stage_1
│           │   ├── genome_scope
│           │   ├── mer_counts.csv
│           │   ├── mer_counts.gs.log.png
│           │   ├── mer_counts.gs.png
│           │   ├── mer_counts.model.txt
│           │   ├── mer_counts.progress.txt
│           │   ├── mer_counts.summary.txt
│           │   ├── out.####.fasta.gz
│           │   ├── out.####.fasta.gz
│           │   ├── ...
│           │   ├── preads_rep.tar.gz
│           │   ├── preads_tan.tar.gz
│           │   ├── raw_reads.###.las
│           │   ├── raw_reads.###.las
│           │   ├── ...
│           │   ├── raw_reads_rep.tar.gz
│           │   ├── raw_reads_tan.tar.gz
│           │   ├── read_counts.csv
│           │   └── read_length_distribution.pdf
│           ├── stage_2
│           │   ├── preads.###.las
│           │   ├── preads.###.las
│           │   └── ...
│           ├── stage_3
│           │   ├── a_ctg.fasta.gz
│           │   ├── a_ctg_tiling_path.gz
│           │   ├── asm.gfa.gz
│           │   ├── asm.gfa2.gz
│           │   ├── contig.gfa2.gz
│           │   ├── ctg_paths.gz
│           │   ├── p_ctg.fasta.gz
│           │   ├── p_ctg_tiling_path.gz
│           │   ├── preads4falcon.fasta.gz
│           │   ├── read_counts.csv
│           │   ├── read_length_distribution.pdf
│           │   ├── sg.gfa.gz
│           │   ├── sg.gfa2.gz
│           │   ├── sg_edges_list.gz
│           │   └── utg_data.gz
│           ├── unzip_stage_1
│           │   ├── ctg_list
│           │   ├── pread_ids.gz
│           │   ├── pread_to_contigs.gz
│           │   ├── rawread_ids.gz
│           │   ├── rawread_to_contigs.gz
│           │   ├── read_to_contig_map.gz
│           │   ├── reads_fastas_##.tar.gz
│           │   ├── reads_fastas_##.tar.gz
│           │   ├── ...
│           │   ├── ref_fastas_###.tar.gz
│           │   ├── ref_fastas_###.tar.gz
│           │   └── ...
│           ├── unzip_stage_2
│           │   ├── all_phased_reads.gz
│           │   ├── rid_to_phase.###.tar.gz
│           │   ├── rid_to_phase.###.tar.gz
│           │   └── ...
│           ├── unzip_stage_3
│           │   ├── all_h_ctg.fa.gz
│           │   ├── all_h_ctg_edges.gz
│           │   ├── all_h_ctg_ids.gz
│           │   ├── all_p_ctg.fa.gz
│           │   ├── all_p_ctg_edges.gz
│           │   ├── asm.gfa.gz
│           │   └── sg.gfa.gz
│           ├── unzip_stage_4
│           │   ├── contig_groups.json
│           │   ├── mappedread_###.tar.gz
│           │   ├── mappedread_###.tar.gz
│           │   └── ...
│           └── unzip_stage_5
│               ├── cns_h_ctg.fasta.gz
│               ├── cns_h_ctg.fastq.gz
│               ├── cns_p_ctg.fasta.gz
│               └── cns_p_ctg.fastq.gz
└── genomic_data
    └── ...
```

The `unzip_stage_5` folder contains the **Primary contigs** (`cns_p_ctg.fasta.gz`) and the **Alternate haplotigs**  (`cns_h_ctg.fasta.gz`), which will be used in the following steps of the assembly pipeline.

These two files must be renamed using the convention of the VGP pipeline. To do this, first select the respective file and then click the pencil button that appears at the right to edit the name. In this example the name should be `fArcCen1_c1.fasta.gz` for the **Primary contigs** file, and `fArcCen1_c2.fasta.gz` for the **Alternate haplotigs** file.

Next, it is required that every intermediate assembly produced during the pipeline is placed in a specific folder. Move the **c1** and **c2** files to the `intermediates` folder by "drag and drop".

<br/>

**!)** At this point of the pipeline it is important to run several assembly metrics to check that all is going well so far, a practice that is repeated after finalizing each step (see the green dots in the [DNAnexus workflow diagram](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#farccen1-assembly-tutorial) at the beginning of the tutorial). To do this, click the green button `+ Add Data` in your working project and search and select **VGP Tools** in the "Other Project" tab. Search and select the latest version of the applets named **asm_stats** and **busco**, and the workflow named **Evaluation KAT Plot**.

To obtain the standard assembly statistics run the **asm_stats** applet using as input the respective assembly to be evaluated. In addition, click the gear icon and complete the "Genome size (bp)" field with the size of the genome in base pairs, and `fArcCen1` in the "species code" field. Inside the `assembly_vgp_standard_1.6` folder, create a new folder with the name `evaluation`. Create a folder inside `evaluation` with the name of the assembly stage to be evaluated (for example, `c1`) and select it as the output folder. Finally, click `Run as Analysis...` to launch the applet.
You should check for an improvement in the assembly metrics with the progress of the pipeline.

To obtain a measure of the completeness of the assembly it is necessary to run **busco** using as input the respective assembly to be evaluated. In addition, click the gear icon and complete the "Augustus species search" filed with the closest species available to your working species, in this example `zebrafish`. Finally, under `Workflow Actions`, select `Set Output Folder`. Create a new folder with the name `busco` inside the `evalutaion/c1` folder and select it as the output folder for the **busco** applet.
You should check for an improvement in the metrics with the progress of the pipeline.

To run the **Evaluation KAT Plot** workflow, select the `_R1_` and `_R2_` files in the `10x` folder as input for the `Remove gembarcodes from 10x reads` stage, and the **c1** and **c2** assemblies as input for the `File Concatenator` stage. For setting the output, create a folder named `KAT_c1c2` inside the `evaluation` folder. Finally, click `Run as Analysis...` to launch the workflow.
You will compare the obtained plot with subsequent steps of the pipeline.


At this stage, the folder structure should look like this:

```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── evaluation
│   │   ├── c1
│   │   │   ├── ...
│   │   │   └── busco
│   │   │       └── ...
│   │   ├── c2
│   │   │   └── ...
│   │   └── KAT_c1c2
│   │       └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       └── fArcCen1_c2.fasta.gz
└── genomic_data
    └── ...
```

**Transferring to S3:** After being sure that each step finished correctly, the stats were checked and the files placed in their respective correct folders, it is a good practice to move the data to the VGP storage in AWS. The data will transfer and a symbolic link will be created to keep files functional and accesible. 
In your working project, click the menu "TOOLS" and select "Tool Library", next search and select the applet **DNAnexus to VGP S3 Exporter**. Select the files generated in the finished step in order to transfer them.

<br/>

### 2. Purge Dups

Briefly, Purge Dups identifies heterozygous regions (mistakenly present in the primary contigs instead of the alternate), and removes them from the primary contigs to place them with the alternate contigs.

In your working project, click the green button `+ Add Data` and search and select **VGP Tools** in the "Other Project" tab. Search and select the latest version of the **Scaffold 1 purge_dups** workflow and click the green button `Add Data`, after which a dialogue box will pop up with a progress bar indicating that the workflow has been copied to the current location of your working project.

Click the workflow to open it in _Run_ mode and create an editable copy of it in case anything is misconfigured or needs to be re-run.

![Purge Dups workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/purge_dups_workflow.png)

The inputs for the workflow are:
* For the `purge_dups` stage: the **c1** file `fArcCen1_c1.fasta.gz` and the _PacBio Reads converted to fasta_ from the `bam_to_fasta` folder
* For the `concat c2+p2` stage: the **c2** file `fArcCen1_c2.fasta.gz`

Finally, under `Workflow Actions`, select `Set Output Folder`. Select `intermediates` as the output folder for the **Scaffold 1 purge_dups** workflow.

All stages of the workflow should now be in the "Runnable" state. Before running, make sure to save your workflow changes (including input specification) by selecting `Workflow Actions` and selecting `Update workflow with changes`. This will make it easier to modify and relaunch the workflow should any failures occur. Finally, click `Run as Analysis...` to launch the workflow.

Once finished, the final output should look like this:
```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── evaluation
│   │   └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       ├── fArcCen1_c2.fasta.gz
│       └── purge_dups
│           ├── haplotig_purged
│           │   ├── fArcCen1_p2.fasta.gz
│           │   ├── fArcCen1_q2.fasta.gz
│           │   └── ...
│           ├── primary_purged
│           │   ├── fArcCen1_p1.fasta.gz
│           │   └── ...
│           ├── read_counts.csv
│           ├── ...
│           ├── read_length_distribution.pdf
│           └── ...
└── genomic_data
    └── ...
```

The **Purged primary** contigs should be contained in the file `fArcCen1_p1.fasta.gz`, and the **Alternate combined** haplotigs should be contained in the `fArcCen1_q2.fasta.gz` file.
Remember to move the `p1` and `q2` files to the `intermediates` folder by "drag and drop".

**!)** Remember to run the required assembly metrics for this stage. You should see an improvement that reflects the performance of the haplotigs purging step when comparing the **KAT** plot obtained for _c1c2_ and for _p1q2_. In addition, check for differences in the **busco** metrics between `c1` and `p1`.

Note: if the _Falcon and Unzip_ step was already run and the **c1** and **c2** are present in the `intermediates` folder but the `bam_to_fasta` folder is not present, you should run the applet **PacBio BAM to FASTA** which can be found by clicking the green button `Start Analysis`. The input of this applet are the _PacBio Sequel Reads_ from the `pacbio` folder and you should set up an output folder named `bam_to_fasta` inside the `intermediates` folder.

<br/>

## Scaffolding

### 1. 10X

The next step of the pipeline consist in two rounds of scaffolding using the 10X Genomics raw reads. To start, click the green button `+ Add Data` in your working project and search and select **VGP Tools** in the "Other Project" tab. Search and select the latest version of the **scaffold_2_scaff10X** workflow and click the green button `Add Data`, after which a dialogue box will pop up with a progress bar indicating that the workflow has been copied to the current location of your working project.

The inputs for the workflow are:
* `assemble_genome_fastagz`: the _Purged primary_ contigs file **p1** from Purge Dups (`fArcCen1_p1.fasta.gz`)
* `scaff_R1_fastqgz`: all reads containing `_R1_*.fastq` in the `10x` folder
* `scaff_R2_fastqgz`: all reads containing `_R2_*.fastq` in the `10x` folder

In addition, specify the `Output Folder` for the workflow to `intermediates` as before. Launch the analysis.

![Configured Scaff10x workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/Scaff10XRunnable.png)

The final output of scaff10x should look like this:
```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── evaluation
│   │   └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       ├── fArcCen1_c2.fasta.gz
│       ├── fArcCen1_p1.fasta.gz
│       ├── fArcCen1_q2.fasta.gz
│       ├── purge_dups
│       │   └── ...
│       └── scaff10x
│           ├── round1
│           │   ├── read-BC_1.fastq.gz
│           │   ├── read-BC_2.fastq.gz
│           │   └── scaffolds.fasta.gz
│           └── round2
│               ├── read-BC_1.fastq.gz
│               ├── read-BC_2.fastq.gz
│               └── scaffolds.fasta.gz
└── genomic_data
    └── ...
```

The output under `round2/scaffolds.fasta.gz` corresponds to `fArcCen1_s1.fasta.gz` using the VGP naming convention.
Remember to move the `s1` file to the `intermediates` folder by "drag and drop".

<br/>

### 2. Bionano

The next step uses the Bionano assembled CMAP data together with the scaffolds file `fArcCen1_s1.fasta.gz` from the previous step to perform hybrid scaffolding on the primary haplotig.

Before copying the workflow, we need to first take a look at the Bionano input data:
```
fArcCen1
└── genomic_data
    ├── bionano
    │   ├── fArcCen1_Saphyr_DLE-1.cmap.gz
    │   ├── fArcCen1_Saphyr_DLE-1_1265239.bnx.gz
    │   └── fArcCen1_Saphyr_DLE-1_1265240.bnx.gz
    ├── ...
    └── ...  
```

From the input data, we can see that there is a single `*.cmap.gz` file generated using the `DLE-1` enzyme. Therefore, copy the **scaffold_3_bionano_1enzyme** workflow from VGP tools as explained before. In some cases you may see two `*.cmap.gz` files in your input data corresponding to two enzymes used for the Bionano data, and therefore you will need to use the _2 enzyme_ workflow.

![Bionano workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/bionano_workflow.png)

The inputs for the workflow are:
* `Merged refineFinal CMAP file`: the `fArcCen1_Saphyr_DLE-1.cmap.gz` from `genomic_data/bionano/` folder
* `FASTA or CMAP file from NGS`: the scaffolds file `fArcCen1_s1.fasta.gz`

As before, specify the workflow output folder to `intermediates`. A `bionano` folder will be automatically created under that.

The Bionano tool is prone to hanging if the memory requirements are not met. Therefore, make sure it is running on a large memory instance, such as the `mem3_ssd1_x32` instance. The app is configured to aggressively time out if it does not complete in under 6 hours. If this happens, rerun the analysis on a larger memory instance. Remember that you can always ask your doubts in the [Slack channel](https://genomeark.slack.com/archives/CE7FU8YAC).


Once the app completes, the output should look smilar as follows:
```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── evaluation
│   │   └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── bionano
│       │   ├── bn_pre_cut_projected_ngs_coord_annotations.bed
│       │   ├── conflicts.txt
│       │   ├── conflicts_cut_status.txt
│       │   ├── conflicts_cut_status_CTTAAG.txt
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_BNGcontigs_HYBRID_SCAFFOLD.xmap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_BNGcontigs_HYBRID_SCAFFOLD_q.cmap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_BNGcontigs_HYBRID_SCAFFOLD_r.cmap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_HYBRID_SCAFFOLD.cmap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_HYBRID_SCAFFOLD_log.txt
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.agp
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.fasta
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.gap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD.xmap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD_NCBI.fasta
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD_NOT_SCAFFOLDED.fasta
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD_q.cmap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD_r.cmap
│       │   ├── fArcCen1_Saphyr_DLE1_bppAdjust_cmap_scaffolds_fasta_NGScontigs_HYBRID_SCAFFOLD_trimHeadTailGap.coord
│       │   ├── fArcCen1_s2.fasta.gz
│       │   ├── hybrid_scaffold_output.tar.gz
│       │   └── ngs_pre_cut_annotations.bed
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       ├── fArcCen1_c2.fasta.gz
│       ├── fArcCen1_p1.fasta.gz
│       ├── fArcCen1_q2.fasta.gz
│       ├── fArcCen1_s1.fasta.gz
│       ├── purge_dups
│       │   └── ...
│       └── scaff10x
│           └── ...
└── genomic_data
    └── ...
```

Remember to move the **s2** file `fArcCen1_s2.fasta.gz` to the `intermediates` folder by "drag and drop".

**!)** During the scaffolding steps like this, besides an improvement in the assembly metrics, it is expectable a decrease in the _contig N50_ due to contig break down.

<br/>

### 3. Salsa

Salsa scaffolding uses Hi-C data to scaffold the hybrid assembly from Bionano.

Take a look at the HiC input data: It will be located under the `genomic_data` folder with the name of the HiC provider who generated it, such as `phase`, `arima` or `dovetail`.

```
fArcCen1
└── genomic_data
    ├── 10x
    │   └── ...
    ├── bionano
    │   └── ...
    ├── pacbio
    │   └── ...
    └── phase
        ├── fArcCen1_DDHiC_R1.fastq.gz
        ├── fArcCen1_DDHiC_R2.fastq.gz
        ├── fArcCen1_S3HiC_R1.fastq.gz
        ├── fArcCen1_S3HiC_R1.fastq.gz
        └── re_bases.txt
```

In addition to the input files, you will need to know the restriction enzymes used to generate the data. For `fArcCen1`,
the sequences are `GATC` since the restriction enzyme employed was MboI. This infomation is reported in the `re_bases.txt` file in the folder with the HiC reads.  

Before starting with the Salsa step, if more than one pair of reads are present in the HiC folder (like in this example), all the _R1_ reads must be concatenated in one file, while the _R2_ reads must be concatenated in other. To do this, click the green button `Start Analysis` in your working project, and search and select the **File Concatenator** applet. The input of the applet are all the _R1_ files in the `phase` folder. Next, under `Workflow Actions`, select `Set Output Folder`, and specify the `phase` folder as output folder. Click `Run as Analysis...` to launch the applet. Finally, you need to **repeat this proceeding for the _R2_ reads** (please be sure of selecting the _R2_ reads **in the same order** the _R1_ reads were selected).

Copy the latest version of the **scaffold_4_salsa** workflow from VGP tools into your project as explained before. The workflow performs the following steps:

1. Align the _HiC reads_ using the Arima mapping pipeline
2. Run _Salsa2_ on aligned reads and the `s2.fasta.gz` scaffolds file
3. Concatenate the output with the _Haplotigs_ from _FALCON Unzip_

![New Salsa workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/new_salsa_workflow.png)

The inputs for the workflow are:
* For the `BWA FASTA Indexer` stage: the scaffolds file `fArcCen1_s2.fasta.gz` for the `Genome` input
* For the `Arima mapping` stage: the HiC reads from the `phase` folder, (concatenated) _R1_ reads file for the `Forward Reads` input and (concatenated) _R2_ reads file for the `Reverse Reads` input.
* For the `Salsa` stage: select the gear icon and specify the HiC restriction enzyme (`GATC`) as the "Restriction enzyme bases" input
* For the `concat s3+q2+mito` stage: the **Alternate combined** haplotigs contained in the file `fArcCen1_q2.fasta.gz` for `q2 input`. If a mitogenome is available for this species, it should be incorporated as input in this stage.

IMPORTANT: You should check for the mitogenome availability in the [VGP GenomeArk website](https://vgp.github.io/genomeark). If it is present, you need to download to your computer, and upload to your DNAnexus project. To do this, click the green button `+ Add Data` and select the file from your computer. If the mitogenome is not available, please ask in in the ["training" channel of the VGP Slack](https://genomeark.slack.com/archives/CE7FU8YAC).


In addition, click the gear icon next to the File Concatenator app to specify the output name: `fArcCen1_s4.fasta.gz`. Remember to configure `intermediates/salsa` as the output folder and save the workflow copy before launch the analysis.

The final output should look like this:

```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── evaluation
│   │   └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── bionano
│       │   └── ...
│       ├── hic
│       │   ├── fArcCen1_s3.fasta.gz
│       │   ├── fArcCen1_s4.fasta.gz
│       │   ├── ...
│       │   ├── ...
│       │   ├── ...
│       │   ├── ...
│       │   ├── ...
│       │   ├── ...
│       │   └── ...
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       ├── fArcCen1_c2.fasta.gz
│       ├── fArcCen1_p1.fasta.gz
│       ├── fArcCen1_q2.fasta.gz
│       ├── fArcCen1_s1.fasta.gz
│       ├── fArcCen1_s2.fasta.gz
│       ├── purge_dups
│       │   └── ...
│       └── scaff10x
│           └── ...
└── genomic_data
    └── ...
```

Remember to move the **s3** and **s4** files to the `intermediates` folder by "drag and drop".

<br/>

## Polishing

### 1. Arrow

The next step of the pipeline consist in a polishing of the scaffolds using PacBio data. To start, click the green button `+ Add Data` in your working project and search and select **VGP Tools** in the "Other Project" tab. Search and select the latest version of the **Scaffold 5 Arrow Polish** workflow and click the green button `Add Data`, after which a dialogue box will pop up with a progress bar indicating that the workflow has been copied to the current location of your working project.

The workflow performs the following steps:
1. Align the _PacBio reads_ to the assembly using _Minimap2_
2. Polish of scaffolds using _Arrow_

![Arrow workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/arrow_workflow.png)

The inputs for the workflow are:
* `Reads`: the _PacBio Sequel Reads_ (`subreads.bam`) from the `pacbio` folder
* `Reference genome`: the scaffolds file `fArcCen1_s4.fasta.gz`

Remember to configure `intermediates` as the output folder and save the workflow copy before launch the analysis.

The final output should look like this:

```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── evaluation
│   │   └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── arrow
│       │   ├── mapped_reads
│       │   │   └── ...
│       │   └── polished_output
│       │       ├── fArcCen1_t1.fasta.gz
│       │       └── ...
│       ├── bionano
│       │   └── ...
│       ├── hic
│       │   └── ...
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       ├── fArcCen1_c2.fasta.gz
│       ├── fArcCen1_p1.fasta.gz
│       ├── fArcCen1_q2.fasta.gz
│       ├── fArcCen1_s1.fasta.gz
│       ├── fArcCen1_s2.fasta.gz
│       ├── fArcCen1_s3.fasta.gz
│       ├── fArcCen1_s4.fasta.gz
│       ├── purge_dups
│       │   └── ...
│       └── scaff10x
│           └── ...
└── genomic_data
    └── ...
```

**!)** At this point, you  should check for an increase in the _contig N50_ when compared with **s3**, as a sign of the gap-filling efficiency.

Remember to move the **t1** file to the `intermediates` folder by "drag and drop".

<br/>

### 2. Freebayes

Copy the latest version of the **Scaffold 6 Longranger Freebayes Polish** workflow from **VGP tools** into your project as explained before.

![Longranger Freebayes workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/longranger_freebayes_workflow.png)

The inputs for the workflow are:
* For the `10X Longranger Reference Builder` stage: the polished scaffolds file `fArcCen1_t1.fasta.gz`
* For the `10X Longranger Align Only Analysis` stage: all the `fastq` reads in the `10x` folder

Remember to configure `intermediates` as the output folder and save the workflow copy before launch the analysis.


The final output should look like this:

```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── evaluation
│   │   └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── arrow
│       │   └── ..
│       ├── bionano
│       │   └── ...
│       ├── hic
│       │   └── ...
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       ├── fArcCen1_c2.fasta.gz
│       ├── fArcCen1_p1.fasta.gz
│       ├── fArcCen1_q2.fasta.gz
│       ├── fArcCen1_s1.fasta.gz
│       ├── fArcCen1_s2.fasta.gz
│       ├── fArcCen1_s3.fasta.gz
│       ├── fArcCen1_s4.fasta.gz
│       ├── fArcCen1_t1.fasta.gz
│       ├── longranger_freebayes_round_1
│       │   ├── 10x_longranger
│       │   │   └── ...
│       │   ├── freebayes
│       │   │   ├── fArcCen1_t2.fasta.gz
│       │   │   └── ...
│       │   └── QV
│       │       └── ...
│       ├── longranger_freebayes_round_2
│       │   ├── 10x_longranger
│       │   │   └── ...
│       │   ├── freebayes
│       │   │   ├── fArcCen1_t3.fasta.gz
│       │   │   └── ...
│       │   ├── QV
│       │   │   ├── qv_report.txt
│       │   │   └── ...
│       │   └── stats
│       │   │   ├── fArcCen1_alt.asm.20191231.fasta.gz
│       │   │   ├── fArcCen1_pri.asm..20191231.fasta.gz
│       │   │   └── ...
│       ├── purge_dups
│       │   └── ...
│       └── scaff10x
│           └── ...
└── genomic_data
    └── ...
```

**!)** You should check for an increase in the QV value after each round, when comparing the `qv_report.txt` file from the folders `longranger_freebayes_round_1/QV` and `longranger_freebayes_round_2/QV`.

Remember to move the **t2** and **t3** files to the `intermediates` folder by "drag and drop".

Move the **final primary and alternative assembly** files to the `assembly_vgp_standard_1.6` folder by "drag and drop".

<br/>

## Final Checks

The final folder structure should look like this:

```
fArcCen1
├── assembly_vgp_standard_1.6
│   ├── fArcCen1_alt.asm.20191231.fasta.gz
│   ├── fArcCen1_pri.asm.20191231.fasta.gz
│   ├── evaluation
│   │   └── ...
│   ├── meryl_genomescope
│   │   └── ...
│   └── intermediates
│       ├── arrow
│       │   └── ..
│       ├── bionano
│       │   └── ...
│       ├── hic
│       │   └── ...
│       ├── falcon_unzip
│       │   └── ...
│       ├── fArcCen1_c1.fasta.gz
│       ├── fArcCen1_c2.fasta.gz
│       ├── fArcCen1_p1.fasta.gz
│       ├── fArcCen1_q2.fasta.gz
│       ├── fArcCen1_s1.fasta.gz
│       ├── fArcCen1_s2.fasta.gz
│       ├── fArcCen1_s3.fasta.gz
│       ├── fArcCen1_s4.fasta.gz
│       ├── fArcCen1_t1.fasta.gz
│       ├── fArcCen1_t2.fasta.gz
│       ├── fArcCen1_t3.fasta.gz
│       ├── longranger_freebayes_round_1
│       │   └── ...
│       ├── longranger_freebayes_round_2
│       │   └── ...
│       ├── purge_dups
│       │   └── ...
│       └── scaff10x
│           └── ...
└── genomic_data
    └── ...
```

This is a typical example of how contig and scaffold N50 should behave during the assembly process:

![N50 plot](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/contig_scaffold_size_plot.png)

<br/>

Finally, it is requested to generate a _HiC heatmap_, a _KAT plot_, and a _BUSCO search_ in the final assembly before sending the genome to curation. In addition, it is also requested to record relevant pipeline information that is required when the genome assembly is submitted to the NCBI and EBI archives. Luckily, all this tasks can be done in a single workflow!

Copy the latest version of the **Presubmission** workflow from **VGP tools** into your project as explained before. Click the workflow to open it in _Run_ mode and under `Workflow Actions`, select `Set Output Folder`. Create a new folder with the name `Presubmission` inside the `evalutaion` folder and select it as the output folder for the **Presubmission** workflow.

![Presubmission workflow](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/presubmission_workflow.png)

The inputs for the workflow are:

* In the `BWA FASTA Indexer` stage: the **pri.asm** file `fArcCen1_pri.asm.20191231.fasta.gz` for the `Genome` input
* In the `Arima mapping` stage: the HiC reads from the `phase` folder, (concatenated) _R1_ reads file for the `Forward Reads` input and (concatenated) _R2_ reads file for the `Reverse Reads` input (remember that concatenation was necessary if more than one pair of reads were present in the HiC folder, check again the [Salsa step](https://github.com/VGP/vgp-assembly/blob/master/tutorials/DNAnexus_workflow_1.6_tutorial.md#3-salsa) for clarification).

Under the stages `BWA FASTA Indexer`, `Arima mapping`, and `pretext`, click the gear icon to open the parameters panel and fill the "Output Folder" parameter with `pretext`

* For the `Remove gembarcodes from 10x reads` stage: the `_R1_` and `_R2_` files in the `10x` folder
* For the `File Concatenator` stage: the **pri.asm** file `fArcCen1_pri.asm.20191231.fasta.gz`, and the **alt.asm** file `fArcCen1_alt.asm.20191231.fasta.gz`

Under the stages `Remove gembarcodes from 10x read`, `File Concatenator`, and `kat`, click the gear icon to open the parameters panel and fill the "Output Folder" parameter with `kat`

* For the `BUSCO genome assembly quality assessment` stage, click the gear icon and complete the "Augustus species search" filed with the closest species available to your working species, in this example `zebrafish`. In addition, fill the "Output Folder" parameter with `busco`

* For the `sw_version` stage, click the gear icon to open the parameters panel and fill the "Output Folder" parameter with `version`


All the stages of the workflow should now be in the "Runnable" state. Click `Run as Analysis...` to launch the workflow.

<br/>

**!)** Once all is finished, please share your genome stats and plots in the [Slack channel](https://genomeark.slack.com/archives/CE7FU8YAC).

<br/>

If everything went well, and after transferring the generated files to S3 with the **DNAnexus to VGP S3 Exporter** applet, you should send a mail to the curation team following a specific template with links and attachments. For this example the mail should look like this:

* Subject:
```
VGP fArcCen1 ready for curation
```
* Body:
```
```

The required QV value can be obtained from `qv_report.txt` file in the `longranger_freebayes_round_2/QV` folder.


<br/>

<br/>

--

This document is being upgraded by D.N. De Panis during a postdoctoral stay at Camila Mazzoni group ([Berlin Center for Genomics in Biodiversity Research](https://begendiv.de)) with the valuable collaboration of VGP members.
