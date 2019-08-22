These steps describe how to call with freebayes and then use `bcftools consensus` to create a polished fasta file. Requirements are:

* freebayes (>=1.1)
* bcftools (>=1.8)
* samtools


0. Setup:

		bam=/path/to/bam/file
		fasta=/path/to/assembly/fasta
		sample=sampleName

1. Index the assembly fasta file and input BAM/CRAM file if not already indexed:

		samtools faidx $fasta
		samtools index $bam

2. Parallelise into chunks across the genome. This awk command will print out the commands to be run:

		awk '{print "freebayes --bam $bam --region "$1":1-"$2" -f $fasta | bcftools view --no-version -Ob -o "$1":1-"$2".bcf"}' $fasta.fai

3. Make a list of the output files (in the same order as the reference) to be concatenated together when done:

		awk '{print $1":1-"$2".bcf"}' $ref.fai > concat_list.txt

4. When all parallel jobs finished, concatenate and normalise non-REF variants:

		bcftools concat -nf concat_list.txt | bcftools view -Ou -e'type="ref"' | bcftools norm -Ob -f $fasta -o $sample.bcf
		bcftools index $sample.bcf

5. Create the polished fasta with `bcftools consensus`. This will do very basic filtering of the callset (`QUAL>1`) and select only homozygous ALT and heterozygous non-REF sites. At heterozygous non-REF sites (both alleles do not match the reference fasta), the longest allele will be chosen.

		bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f $fasta $sample.bcf > $sample.fasta

	Alternatively, if you have a phased VCF and would like to consistently pick one of the haplotypes at all heterozygous sites, set `-H` to `1` or `2`:

		bcftools consensus -i'QUAL>1' -H1 -f $fasta $sample.bcf > $sample.fasta
