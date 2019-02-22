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
    <genome_id>_S1_L001_I1_001.fastq.gz
    <genome_id>_S1_L001_R1_001.fastq.gz
    <genome_id>_S1_L001_R2_001.fastq.gz
    ```
  * bionano [platform=Irys|Saphyr]
    ```
    <genome_id>_<platform>_<enzyme>[_jobid].bnx.gz
    <genome_id>_<platform>_<enzyme>.cmap.gz
    ```
  * arima
    ```
    <genome_id>_<runID>_R1.fastq.gz
    <genome_id>_<runID>_R2.fastq.gz
    re_bases.txt	(ex. GATC,GANTC)
    ```
  * illumina (Optional)
    ```
    <genome_id>_<runID>_R1.fastq.gz
    <genome_id>_<runID>_R2.fastq.gz
    ```

* assembly_\[pipeline\]_\[ver\] (pipeline: vgp_standard, cambridge, ...)
  * intermediates
    * falcon_unzip   &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; FALCON unzip intermediate files
    * purge_haplotigs	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; purge_haplotigs intermediate files
    * scaff10x	&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Scaff10X intermediate files
    * bionano &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Bionano TGH intermediate files
    * salsa &nbsp;&nbsp;&nbsp;&nbsp; Salsa intermediate files
    * arrow &nbsp;&nbsp; Arrow polishing intermediate files
    * longer_freebayes_round1 &nbsp;&nbsp;&nbsp;&nbsp; Longranger freebayes polishing intermediate files (round1)
    * longer_freebayes_round2 &nbsp;&nbsp;&nbsp;&nbsp;&nbsp; Longranger freebayes polishing intermediate files (round2)
    * ...
    ```
    <genome_id>_c1.fasta	Pacbio FALCON-Unzip assembly primary contigs (haplotype 1)
    <genome_id>_c2.fasta	Pacbio FALCON-Unzip assembly associated haplotigs (haplotype 2)
    <genome_id>_p1.fasta	purge_haplotigs curated primary assembly (taking c1 as input)
    <genome_id>_p2.fasta	purge_haplotigs curated haplotigs (purged out from c1)
    <genome_id>_q2.fasta	c2 + q2 for future polishing
    
    <genome_id>_s1.fasta	2-rounds of scaff10x; scaffolding p1.fasta
    <genome_id>_s2.fasta	Bionano TGH; hybrid scaffold of 2 enzymes over s1.fasta
    <genome_id>_s3.fasta	Salsa scaffolding with Arima hiC libraries over s2.fasta
    <genome_id>_t1.fasta	Arrow polishing over s3.fasta
    <genome_id>_t2.fasta	1 round of longranger_freebayes polishing over t1.fasta
    <genome_id>_t3.fasta	2nd round of longranger_freebayes polishing over t2.fasta
    ```

  ```
  <genome_id>.pri.asm.YYYYMMDD.fasta.gz       Final assembly (primary)
  <genome_id>.alt.asm.YYYYMMDD.fasta.gz       Final assembly (alternate haplotigs)
  ```
* assembly_curated
  ```
  <genome_id>.pri.cur.YYYYMMDD.fasta.gz       Final curated assembly (primary)
  <genome_id>.alt.cur.YYYYMMDD.fasta.gz       Final curated assembly (alternate haplotigs)
  <genome_id>.pri.cur.YYYYMMDD.agp            Chromosome assignments for <genome_id>.pri.cur.YYYYMMDD.fasta.gz
  <genome_id>.pri.cur.YYYYMMDD.MT.fasta.gz    Mitochondrial genome assembly (optional)
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
      <genome_id>_R1.fastq.gz
      <genome_id>_R2.fastq.gz
      ```

## Detailed intermediate assembly names and rules for v1

| intermediate.fasta	| full_verbal | description |
|:------------- | :---------- | :-----------|
|c1.fasta	| pac_fcnz_hap1	| pac_fcnz_hap1: Pacbio FALCON-Unzip assembly primary contigs |
|c2.fasta	| pac_fcnz_hap2	| pac_fcnz_hap2: Pacbio FALCON-Unzip assembly alternate haplotigs |
|p1.fasta	| pac_fcnz_hap1_purg_prim	| prim: purge_haplotigs curated primary |
|p2.fasta	| pac_fcnz_hap1_purg_alt	| purg: purged haplotigs |
|s1.fasta	| pac_fcnz_hap1_10x_scaff10x	|scaff10x: 2-rounds of scaff10x |
|s2.fasta	| pac_fcnz_hap1_10x_scaff10x_bio_tgh	|tgh: bionano TGH; hybrid scaffold of 2 enzymes. *Make sure to include the NOT_SCAFFOLDED leftovers.*|
|s3.fasta	| pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa | arim_salsa: maximum 5-round of Salsa scaffolding from Arima hiC libraries |
|t1.fasta	| pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa_arrow	| arrow: arrow polishing with gap filling |
|t2.fasta |	pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa_arrow_frb1 |	longranger + freebayes polishing round 1 |
|t3.fasta |	pac_fcnz_hap1_10x_scaff10x_bio_tgh_arim_salsa_arrow_frb2 |	longranger + freebayes polishing round 2 |

