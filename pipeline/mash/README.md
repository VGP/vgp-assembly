# Mash screen and dist

## Get distances in your sequence sets

Let’s get distance among your sequence sets.
If any contamination occurred in one particular sequencing platform, it will be visibly clustering themselves.

`sketch.sh` will go through all sequence fasta/bam files and create sketches.

``` 
mash sketch -k 21 -s 10000 -r -m 1 -o $output_msh $input_fasta
```

Combine all the sketches and get the distance.
```
mash paste combined.msh *.msh
mash dist -t combined.msh combined.msh > combined.tbl
```

The plotting script will provide a visual glance on the clustering result.
```
Rscript plot.R $out
```


## Screen what is in your sequence

Download and run with pre-sketched genomes: [refseq.genomes.k21s1000.msh](https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh).

And run mash screen against that reference genome set.

```
mash screen -w refseq.genomes.k21s1000.msh $input_read > out
```

Interpreting the result
The output file contains
[identity, shared-hashes, median-multiplicity, p-value, query-ID, query-comment]

Usually, closely related species (or the same species) will be popped up with identity>85% and high p-value.
So, anything falling in this boundary, but not expected may be a contaminated sequence.

However, there is also a high chance that the genomes randomly share some repetitive sequences. Those may also be popped up in the top.
Also, we have frequently seen a phage sequence used for illumina library preparation being popped up on top.
This is: NC_001422.1 Enterobacteria phage phiX174 sensu lato, complete genome

### For your interest..
For more details about mash screen, see this blog post: [ what’s in my sequencing run?](https://genomeinformatics.github.io/mash-screen/).

For preparing sketches your self, use
```
mash sketch -k 21 -s 10000 -r -m 1 genome.fasta
```


