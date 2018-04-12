# Data structure on DNAnexus and AWS GenomeArk


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
    ```
    <movie>.subreads.bam
    <movie>.subreads.bam.pbi
    <movie>.scraps.bam (Optional for running QC plots)
    ```
  * 10x
    ```
    <sample>_S1_L001_I1.fastq.gz
    <sample>_S1_L001_R1.fastq.gz
    <sample>_S1_L001_R2.fastq.gz
    ```
  * bionano [platform=Irys|Saphyr]
    ```
    <sample>_<platform>_<enzyme>[_jobid].bnx.gz
    <sample>_<platform>_<enzyme>.cmap.gz
    ```
  * arima
    ```
    <sample>_<runID>_R1.fastq.gz
    <sample>_<runID>_R2.fastq.gz
    re_bases.txt	(ex. GATC,GANTC)
    ```
  * illumina (Optional)
    ```
    <sample>_<runID>_R1.fastq.gz
    <sample>_<runID>_R2.fastq.gz
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
      ```
      <movie>.subreads.bam
      <movie>.subreads.bam.pbi
      <movie>.subreadset.xml
      ```
    * illumina
      ```
      <sample>_R1.fastq.gz
      <sample>_R2.fastq.gz
      ```

## Detailed intermediate assembly names and rules for v1

| intermediate.fasta	| full_verbal | description |
|:------------- | :---------- | :-----------|
|c1.fasta	| pac_fcnz_hap1	| pac_fcnz_hap#: Pacbio FALCONunzip assembly hap# |
|c2.fasta	|pac_fcnz_hap2	||
|s1.fasta	|pac_fcnz_hap1_10x_scaff10x	|scaff10x: 2-rounds of scaff10x joining pac_fcnz_hap1 and 10x_spnv_hap1
|s2.fasta	|pac_fcnz_hap1_10x_scaff10x_bio_tgh	|tgh: bionano TGH; hybrid scaffold of 2 enzymes. *Make sure to include the NOT_SCAFFOLDED leftovers.*|
|s3.fasta	|pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa | arim_salsa: 1-round of Salsa scaffolding from Arima hiC libraries |
|t1.fasta	|pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa_bio_tgh_all	|bio_tgh: bio_tgh with a space of manual curation (cutting off severe scaffolding mistakes), and concatinating the pac_fcnz_hap2.|
|t2.fasta|	pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa_bio_tgh_all_pbjl|	pbjelly |
|t3.fasta|	pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa_bio_tgh_all_pbjl_arrow|	arrow |
|t4.fasta	|pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa_bio_tgh_all_pbjl_arrow_pilon2	|2 rounds of pilon with 10X illumina reads |
