# Minimap2 + purge_dups
Running purge_dups with Minimap2

## Requirements
* [Samtools](https://github.com/samtools/samtools) in path
* [Minimap2](https://github.com/lh3/minimap2) in path
* [purge_dups](https://github.com/dfguan/purge_dups) installed as `$tools/purge_dups/bin/purge_dups`

## Workflow
```
Usage: ./_submit_purge_dups.sh <asm.fasta> <input.fofn>
```

### Index
Generate `$ref.idx` so we don’t have to compute it every time the alignment runs with `minimap2_idx.sh`
```
minimap2_idx.sh asm.fasta
```

### Align
Generate `.paf.gz` files for each `read.fastq.gz` in `input.fofn` with `minimap2.sh`
```
minimap2 -x $preset -t $cpus $ref.idx $qry | gzip -c - > $out.read.paf.gz
```

Generate self-alignment `self.paf.gz` of the `asm.fasta` with `minimap2_self.sh`
```
split_fa $asm > $asm.split
minimap2 -xasm5 -DP -t $cpus $asm.split $asm.split | gzip -c - > $asm.split.self.paf.gz
```

### Purge
Gather the `.paf.gz` files, calculate `cutoffs`, and purge duplications
```
pbcstat *.read.paf.gz
calcuts PB.stat > cutoffs 2>calcults.log
purge_dups -2 -T cutoffs -c PB.base.cov $asm.split.self.paf.gz > dups.bed 2> purge_dups.log
get_seqs dups.bed $asm > purged.fa 2> hap.fa 
```

## Output
* purged.fa : assembly to use for the next scaffolding step
* hap.fa : purged out alternate haplotype and junk sequences

## Notes
Check the `cutoffs` and `PB.stat` to make sure the cutoffs are properly chosen.
Run the last `purge_dups` with modified cutoffs when necessary.
