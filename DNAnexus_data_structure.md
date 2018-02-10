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
        subreads.bam
        subreads.bam.pbi
        ```
    * bam_scraps
        ```
        scraps.bam
        ```
    * fasta (or fastq)
        ```
        subreads.fasta
        ```
  * 10X
    * fastq
      ```
      <sample>_S1_L001_I1_001.fastq
      <sample>_S1_L001_R1_001.fastq
      <sample>_S1_L001_R2_001.fastq
      ```
  * bionano
    * bnx
      ```
      bspqi.bnx
      bsssi.bnx
      ```
    * cmap
      ```
      bspqi.cmap
      bsssi.cmap
      ```
  * arima
    *	fastq
         ```
         R1.fastq
         R2.fastq
         re_bases.txt	(ex. GATC,GANTC)
         ```
   * illumina (Optional)
     * fastq
       ```
       R1.fastq
       R2.fastq
       ```

* assembly_\[ver\]
  * intermediates
    * fcnz   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; FALCON unzip intermediate files
    * scaff10x	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Scaff10X intermediate files
    * tgh &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Bionano TGH intermediate files
    * salsa &nbsp;&nbsp;&nbsp;&nbsp; Salsa intermediate files
    * pbjelly &nbsp;&nbsp; PBJelly intermediate files
    * longr &nbsp;&nbsp;&nbsp;&nbsp; Longranger intermediate files
    * freebayes &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Pilon intermediate files
    * ...
    ```
    <sample>_<ver>_c1.fasta	Pacbio FALCON-Unzip assembly primary contigs (haplotype 1)
    <sample>_<ver>_c2.fasta	Pacbio FALCON-Unzip assembly associated haplotigs (haplotype 2)
    <sample>_<ver>_s1.fasta	2-rounds of scaff10x; scaffolding c1.fasta
    <sample>_<ver>_s2.fasta	Bionano TGH; hybrid scaffold of 2 enzymes over s1.fasta
    <sample>_<ver>_s3.fasta	1-round of Salsa scaffolding with Arima hiC libraries over s2.fasta
    <sample>_<ver>_t1.fasta	Bionano TGH + manual curation over s3.fasta
    <sample>_<ver>_t2.fasta	PBJelly over t1.fasta
    <sample>_<ver>_t3.fasta	polishing with 10X reads over t2.fasta
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
