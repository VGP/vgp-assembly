GenomeScope plots overview
==========================
Based on [Interpreting GenomeScope profiles for VGP genome assemblies](https://hackmd.io/@tlama/Sk1HmluTH) by Tanya Lama.

[GenomeScope](https://github.com/schatzlab/genomescope) is used to estimate the overall characteristics of a genome including genome size, heterozygosity rate and repeat content from unprocessed short reads using a kmer-based statistical approach.


The results are summarized in two plots, where each plot is kmer coverage (x) by kmer frequency (y). Therefore, a point in the plot represents how many different kmers spans a given coverage. As an example, repetitive regions are usually composed of particular kmers, so there are not going to be many different kmers with a high coverage (i.e. one single kmer is highly repeated and therefore has high coverage). In the VGP pipeline, a kmer of 31bp long is used for the meryl+genomescope workflow.


The main result is the genome size estimation, which in this example can be abbreviated in 1,273 Gbp. Additionally, the interpretation of the profiles allows identifying where haploid and diploid peaks are expected, but also to have an idea of the sequencing performance. Shortly, the average kmer coverage is going to be similar to the effective coverage of the reads in the genome, which is useful as an estimation of the sequencing depth.
![Plot genomescope 1](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/genomescope_plot.png)


The log plot allows to represent those 31bp long kmers that are highly repeated and therefore have high coverage (e.g. repetitive regions). 
![Plot genomescope 2](https://github.com/VGP/vgp-assembly/blob/master/tutorials/images_1.6/genomescope_log_plot.png)

<br/>

Remember that you can always comment your results and ask your doubts in the ["training" channel of Slack](https://genomeark.slack.com/archives/CE7FU8YAC)
