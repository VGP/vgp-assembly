<!-- dx-header -->
# Scaff10x: The genome assembly pipeline based on read clustering (DNAnexus Platform App)

The genome assembly pipeline based on read clustering

## How should I set the parameter?

a. SMALT is notably slower than BWA. So your first try goes to BWA.

b. The block value is very important. The default value of 2500 is very conservative and you may increase this value to say 5000 or 10000 to improve the length of scaffolds.

c. The default numbers of "reads" and "link" are based on 30X coverage of the data. You need to increase these values if you have say 60X reads.

d. To get the best results, you may run the pipeline twice, here the new genome-assembly.fasta input file as the scaffolded file from last time.

e. By using the option of "-longread 1", the pipeline performs aggressively mapping score filtering on small PacBio/ONT contigs.  

## How does this app work?

This app uses the Scaff10. For general information, consult the manual at:

https://sourceforge.net/projects/phusion2/files/scaff10x/

This app performs the following steps:

- Extract barcode information from fastq data
- Map with BWA or SMALT. The program can use the supplied SAM file as well
- Scaffolding 
- Run Break10X to indentify misassembly

