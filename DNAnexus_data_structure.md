# Assembly data structure on DNAnexus


## Guidance
* Each bulletin points are folder names,
  * with each daughter folders/files shown below

```
File names with extensions are provided in this block.
```

Once we visit a VGP project, the top level folders are expected to have:
* genomic_data
* assembly_\[ver\]
* transcriptomic_data

## Full sub-folders and files

* genomic_data
  * pacbio
    * bam_subreads
        ```
        <movie>.subreads.bam
        <movie>.subreads.bam.pbi
        ```
    * bam_scraps
        ```
        <movie>.scraps.bam
        ```
    * fasta (or fastq)
        ```
        <movie>.subreads.fasta.gz
        ```
  * 10X
    * fastq
      ```
      <sample>_S1_L001_I1_001.fastq.gz
      <sample>_S1_L001_R1_001.fastq.gz
      <sample>_S1_L001_R2_001.fastq.gz
      ```
  * bionano [platform=Irys|Saphyr]
    * bnx
      ```
      <sample>_<platform>_<enzyme>.bnx.gz
      ```
    * cmap
      ```
      <sample>_<platform>_<enzyme>.cmap.gz
      ```
  * arima
    *	fastq
         ```
         <sample>.<runID>.R1.fastq.gz
         <sample>.<runID>.R2.fastq.gz
         re_bases.txt	(ex. GATC,GANTC)
         ```
  * illumina (Optional)
    * fastq
      ```
      <sample>.<runID>.R1.fastq.gz
      <sample>.<runID>.R2.fastq.gz
      ```

* assembly_\[ver\]
  * intermediates
    * fcnz   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; FALCON unzip intermediate files
    * spnv	&nbsp;&nbsp;&nbsp;&nbsp; Supernova intermediate files
    * qm	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Quickmerge intermediate files
    * tgh &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Bionano TGH intermediate files
    * salsa &nbsp;&nbsp;&nbsp;&nbsp; Salsa intermediate files
    * pbjelly &nbsp;&nbsp; PBJelly intermediate files
    * arrow &nbsp;&nbsp;&nbsp; Arrow intermediate files
    * longr &nbsp;&nbsp;&nbsp;&nbsp; Longranger intermediate files
    * pilon &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Pilon intermediate files
    * ...
    ```
    c_f1.fasta	Pacbio FALCONunzip assembly haplotype 1
    c_f2.fasta	Pacbio FALCONunzip assembly haplotype 2
    c_t1.fasta	10X Genomics Supernova assembly haplotype 1
    c_t2.fasta	10X Genomics Supernova assembly haplotype 2
    s_1a.fasta	2-rounds of quickmerge; joining c_f1.fasta and c_t1.fasta
    s_1b.fasta	Bionano TGH; hybrid scaffold of 2 enzymes over s_1a.fasta
    s_1c.fasta	1-round of Salsa scaffolding from Arima hiC libraries over s_1b.fasta
    s_3a.fasta	Bionano tgh over s_1c.fasta + c_f2.fasta
    s_3b.fasta	PBJelly over s_3a.fasta
    s_3c.fasta	Arrow over s_3b.fasta
    s_3d.fasta	2 rounds of pilon with 10X reads over s_3c.fasta
    ```

  ```
  <sample>_<ver>.fasta       Final <ver> (i.e. v1) assembly
  <sample>_<ver>_hap1.fasta				
  <sample>_<ver>_hap2.fasta				
  <sample>_<ver>.agp          optional
  <sample>_<ver>.gfa.bed      pe bed format file using hap1 as 1st pair, hap2 as 2nd pair
  <sample>_<ver>.gfa
  ```

* transcriptomic_data
  * \[tissue\]{#}
    * pacbio
      * bam
        ```
        subreads.bam
        subreads.bam.pbi
        subreadset.xml
        ```
      * illumina
        * fastq
          ```
          R1.fastq
          R2.fastq
          ```

## Detailed intermediate assembly names and rules for v1b

| intermediate.fasta	| full_verbal | description |
|:------------- | :---------- | :-----------|
|c_f1.fasta	| pac_fcnz_hap1	| pac_fcnz_hap#: Pacbio FALCONunzip assembly hap# |
|c_f2.fasta	|pac_fcnz_hap2	||
|c_t1.fasta	|10x_spnv_hap1	|10x_spnv_hap#: 10X Genomics Supernova assembly hap# |
|c_t2.fasta	|10x_spnv_hap2	||
|s_1a.fasta	|pac_fcnz_10x_spnv_qm	|qm: 2-rounds of quickmerge joining pac_fcnz_hap1 and 10x_spnv_hap1
|s_1b.fasta	|pac_fcnz_10x_spnv_qm_bio_tgh	|bio_tgh: bionano TGH; hybrid scaffold of 2 enzymes. *Make sure to include the NOT_SCAFFOLDED leftovers.*|
|s_1c.fasta	|pac_fcnz_10x_spnv_qm_bio_tgh_arim_salsa | arim_salsa: 1-round of Salsa scaffolding from Arima hiC libraries |
|s_3a.fasta	|pac_fcnz_10x_spnv_qm_bio_tgh_arim_salsa_bio_tgh_all	|bio_tgh_all: bio_tgh, this time including pac_fcnz_hap2. Again, don't forget to include the *NOT_SCAFFOLDED leftovers.* |
|s_3b.fasta|	pac_fcnz_10x_spnv_qm_bio_tgh_arim_salsa_bio_tgh_all_pbjl|	pbjelly |
|s_3c.fasta|	pac_fcnz_10x_spnv_qm_bio_tgh_arim_salsa_bio_tgh_all_pbjl_arrow|	arrow |
|s_3d.fasta	|pac_fcnz_10x_spnv_qm_bio_tgh_arim_salsa_bio_tgh_all_pbjl_arrow_pilon2	|2 rounds of pilon with 10X illumina reads |
