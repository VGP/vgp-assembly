# Asset: Reliable blocks

## Reliable blocks
Reliable blocks are regions with support from CLR subreads, linked reads, raw molecule and label information (bnx) of the optical maps, and Hi-C maps from the same individual collected. 

Low or high coverage regions are excluded, which are indicators of mis-assemblies. To overcome mapping biases, this script requires at least two independent platforms to agree for being structurally ‘reliable’. For example, a repeat region longer than a CLR read may cause abnormal high coverage in CLR and linked reads from mapping biases, even if the region was well assembled locally. Longer range data such as optical maps and Hi-C interactions can complement this bias and indicate structural reliance.

Read mapping is performed individually for each platform, and coverage support information is collected on `$ref.fasta` with [Asset v1.0.2](https://github.com/dfguan/asset). Here we include a brief description of each parameter for each platform as well as codes used to generate supporting regions. A manuscript for Asset will follow with more detailed information.

Regions with not enough support (hereby “low support”) were merged when less than 100 bp apart. Cutoff for defining low support differs per platform, as noted below per platforms. Reliable regions were calculated by excluding these low support regions from the assembly. Because sequencing coverage naturally drops at the end of scaffolds for optical maps and Hi-C, we included any low support region as “reliable” that overlaps 1 kbp of each ends in the scaffolds. All available optical maps can be used, including those not used for hybrid scaffolding.


## Requirements
* [Asset](https://github.com/dfguan/asset): Tested with v1.0.2
* [Minimap2](https://github.com/lh3/minimap2): For long-reads, CLR / ONT / HiFi, tested with v2.17
* [BWA](https://github.com/lh3/bwa): For short-reads (10XG or Hi-C)
* [Bionano Solve](https://bionanogenomics.com/support/software-downloads/): For optical maps, tested with vSolve3.3_10252018
* [Samtools](https://github.com/samtools/samtools)
* [BEDTools](https://github.com/arq5x/bedtools2): Tested with v2.92.2

## Setting up

Download [asset](https://github.com/dfguan/asset) and install following the instruction.

Set `$asset` path variable pointing to asset top level directory
```
export asset=/path/to/asset
```

## Workflow
This pipeline does not have a central `_submit_` script, instead is providing submitters for each platform separately. The workflow is as following:

1. Gaps
2. Supportive regions for each platform
3. Reliable blocks

## Input
* Long-reads or alignment (.bam), works for CLR / HiFi / ONT
* 10XG Linked reads or alignment (.bam), accepts bam file generated with [bwa](https://github.com/VGP/vgp-assembly/tree/master/pipeline/bwa) or [longranger](https://github.com/VGP/vgp-assembly/tree/master/pipeline/longranger)
* Illumina WGS reads or alignment
* Bionano .bnx and .cmaps
* Hi-C reads or alignment (.bam), consider using [this pipeline](https://github.com/VGP/vgp-assembly/tree/master/pipeline/pretext)
* Pre-calculated mean coverage for each long-read platform

## 1. Gaps
```
Usage: ./asset.sh <asm.fasta> <name>
    <asm.fasta>: absolute path.
    <name>: <asm.fasta> will be linked as <name.fasta>
```

This step generates 2 important files:
* `gaps.bed` : Gaps
* `asm.ends.bed` : start, end of each scaffold 1 kbp end coordinates

Make directories for each platform and link these before launching asset. 

The above code is running:
```
# Get gaps
$asset/bin/detgaps $name.fasta > gaps.bed
 
# Get scaffold length
samtools faidx $name.fasta
 
# Get scaffold region
awk '{print $1"\t0\t"$2}' $name.fasta.fai > asm.bed
 
# Get scaffold ends while ignoring scaffolds <2kb
cat asm.bed | \
 awk '($3-$2) > 2000 {print $1"\t0\t1000\n"$1"\t"($3-1000)"\t"$3}' \
 > asm.ends.bed
```

## 2. Supportive regions for each platform

### PacBio CLR reads (also applicable for PacBio HiFi or ONT reads)
```
minimap2/_submit_minimap2.sh  <ref.fasta> [genome_cov]
	Assumes <input.fofn> in the same path
	[genome_cov]: set ast_pb -M to genome_cov*2.5
```

Reads are aligned using `minimap2` with `-x map-pb`. Recent ONT base call (Guppy 3.6.0+) showed a similar or lower error rate compared to CLR, which makes it suitable to use the same `-x map-pb` option.
 
Once all `.paf` files are collected, ast_pb runs with `-M $max`, which the maximum threshold identified being from `(sequencing mean coverage) x 2.5`. The mean coverage could be inferred from the estimated haploid genome size / total bases or from the alignments.
 
By default, `ast_pb` only includes read alignments with a minimum of 600 bases. Where `r` is a read, `s` the starting and `e` the ending coordinate of the alignment of `r`, any read alignment with `r(s+300, e-300)` is used to avoid errors at read ends. Regions with a minimum of 10 read alignments are excluded with the default `-m` option.
 
Brief help message is as follows: 
```
Usage: ast_pb [options] <PAF_FILE> ...
Options:
 -m INT minimum coverage [10]
 -M INT maximum coverage [400]
 -l INT bases clipped at start and end coordinates of an alignment [300]
 -h  help
```

The submission code is running `minimap2/minimap2_idx.sh` `minimap2/minimap2.sh` for each input `.fastq.gz` and `minimap2/ast.sh` at the end:
```
# Reference indexing
minimap2 -t $cpus -x map-pb -d $ref.idx $ref.fasta

# Align each qry subread .fasta file to the reference index
minimap2 -x map-pb -t $cpus $ref.idx $qry > $out.paf
 
# Accumulate coverage and exclude low and high coverage
pafs=`ls *.paf`
max=`echo $mean_cov | awk '{printf "%.0f\n", $1*2.5}'`
$asset/bin/ast_pb -M $max $pafs > pb_M.bed"
```

### 10X Genomics Linked reads
```
./_submit_10x.sh <genome> [mean_cov]
Asset will run on aligned.bam from longranger and gaps.bed in this dir.
```

The `aligned.bam` file from `longranger align` are assumed to be in place. To get the $max threshold, `ast_10x` was run in two rounds. The first round was run to get the average molecule coverage. The second round was run with `-C $max`, which is `(average molecule coverage) x 3.5`. Note this is set higher than what was used in CLR coverage as the linked reads were aligned to both haplotypes, thus here the average molecule coverage is closer to the haploid coverage.
 
By default, `ast_10x` requires regions to have at least `0.15 x (average molecule coverage)` or 10 molecules, whichever is higher. This script increases this threshold to use `3.5 x (average molecule coverage)` to avoid overcalling due to coverage fluctuations. Molecules are only considered when the average mapping quality of the reads in it is over 20, with the inferred molecule size being longer than 1 kbp. A molecule requires shared barcodes among at least 20 reads, where any two adjacent reads are less than 20 kbp apart. The maximum number of reads in a barcode is restricted to at most 1 million; however, based on the number of reads in barcodes, most barcodes meet this filtering criteria.
 
Brief help message is as follows:
Usage: ast_10x [options] <GAP_BED> <BAM_FILEs> ...
Options:
  -x BOOL use longranger bam [False]
  -b INT  minimum number of reads for each barcode [20]
  -B INT  maximum number of reads for each barcode [1M]
  -c INT minimum molecule coverage. This or -r will be used, whichever is higher. [10]
  -r FLOAT minimum coverage ratio to the average coverage [.15]
  -C INT maximum coverage [inf]
  -q INT minimum average read mapping quality for each molecule [20]
  -l INT minimum length for a molecule [1000]
 -S INT maximum distance allowed between two adjacent reads with identical barcode to be grouped as a molecule [20000]
  -a INT minimum number of barcodes for each molecule [5]
  -h 	help
 
The submission code is running `longranger/ast.sh`:
```
# First round: accumulate molecule coverage to get the mean
$asset/bin/ast_10x -x gaps.bed aligned.bam > 10x.bed
 
# Avg. molecule coverage and max cutoff
mean_cov=`awk '{sum+=$1*$2; total+=$2} END {printf "%.0f\n", sum/total}' TX.stat`
max=`echo $mean_cov | awk '{printf "%.0f\n", $1*3.5}'`
 
# Second round
$asset/bin/ast_10x -x -C $cutoff $gaps aligned.bam > 10x_C.bed
```

### Bionano Optical map raw molecule (bnx)
```
_submit_bionano.sh <ref.fasta>
Assumes input.fofn in the path.
```
* Preparing .bnx files
Merge the `.bnx` molecules prior to alignment when multiple bnx are available from the same sequencing platform (Irys or Saphyr) and label to match `<genome>_<platform>_<label>.bnx`.

Use the [merge_bnx.sh](https://github.com/VGP/vgp-assembly/blob/master/pipeline/bionano/merge_bnx.sh) script to merge the `.bnx` files.
```
merge_bnx.sh <prefix>
    <prefix>: ex. mLynCan4_Saphyr_BspQI for mLynCan4_Saphyr_BspQI_748455.bnx and mLynCan4_Saphyr_BspQI_792415.bnx
    <output>: <prefix>.bnx
```

The `$ref.fasta` first converted to in-silico reference cmaps for each label (BspQI, BssSI, or DLE1) to align available bnx maps accordingly. 
 
The script aligns each `.bnx` molecule file using `RefAligner in Solve v3.3_10252018` using non-haplotype option of the sequencing platform (Irys or Saphyr), as we align bnx from both haplotypes to a pseudo-haplotype assembly. Molecule coverage was obtained with `ast_bion_bnx` using default options, which requires regions to have at least 10 molecule coverage or `0.5 x (average molecule coverage)`, whichever is higher. When multiple bnx files were used, all molecule coverage is gathered using the `union` function of Asset.
 
Brief help message is as follows:
```
Usage: ast_bion_bnx [options] <REF_CMAP> <QUERY_CMAP> <XMAP> <KEY_FN>
Options:
  -m INT   minimum molecule coverage [10]
  -M INT   maximum molecule coverage [inf]
  -r INT   minimum coverage ratio to mean coverage [.5]
  -s FLOAT minimum alignment confidence [0.0]
  -O STR   output directory [.]
  -h   	help
```
 
The submission code is running `bionano/map_bnx.sh` for each input `.bnx` and `bionano/ast.sh` at the end:
```
# Convert primary reference assembly fasta to cmap
perl $solve_dir/HybridScaffold/10252018/scripts/fa2cmap_multi_color.pl -e $enzyme 1 -i $ref -o $output_dir/fa2cmap
 
# Merge if multiple bnx are available from the same platform and label, prefix is for example mLynCan4_Saphyr_BspQI
$tools/bionano/Solve3.3_10252018/RefAligner/7915.7989rel/RefAligner -if $prefix.list -merge -o $prefix -bnx -stdout -stderr
 
# Align bnx to the reference cmap
python $solve_dir/Pipeline/10252018/align_bnx_to_cmap.py --prefix $enzyme --mol $query_map --ref $ref_cmap --ra $solve_dir/RefAligner/7915.7989rel/ --nthreads $cpus --output $output_dir/align --optArgs $solve_dir/RefAligner/7915.7989rel/optArguments_nonhaplotype_"$platform".xml --pipeline $solve_dir/Pipeline/10252018/
 
# Convert to support regions of this .bnx.
# $rmap_fn, $qmap_fn, $xmap_fn, and $key_fn are the output files
# of the above fa2cmap
$asset/bin/ast_bion_bnx $rmap_fn $qmap_fn $xmap_fn $key_fn \
  > $output_dir/bionano_"$tech"_"$enzyme".bed \
  2>ast_bion_bnx_"$tech"_"$enzyme".log
 
# Merge support regions when multiple enzymes were used
$asset/bin/union bnx_*/bionano_*.bed > bn.bed
```

### Hi-C interaction
Running asset with an existing `.bam` file, ex. from  [this process](https://github.com/VGP/vgp-assembly/tree/master/pipeline/pretext):
```
_submit_hic.sh <name>
	Run with an existing .bam file
  Requires gaps.bed
```

Alternatively, generate alignments with `bwa mem` and run asset:
```
_submit_bwa.sh <ref.fasta> [jid]
	<ref.fasta>: no gz!
	Assumes <input.fofn> and gaps.bed in the same path
```

Coverage information is obtained using ast_hic with default options, which excludes regions with less than seven interactions. An interaction is inferred from the distance of a read pair, using the starting coordinates of each read while excluding N-base gaps. Only interactions less than 15 kbp were considered in coverage to avoid noisy long-range interactions for inferring structural reliability.

Brief help message is as follows:
```
Usage: aa_hic [options] <GAP_BED> <BAM_FILEs>
Options:
  -c INT minimum coverage [7]
  -C INT maximum coverage [inf]
  -q INT minimum alignment quality [0]
  -L INT maximum insertion length, gap excluded [15000]
  -h     help
```

The submission code is running `hic/ast.sh`:
```
# Convert alignments to support information
$asset/bin/ast_hic gaps.bed *.bam > hic.bed
```

## 3. Reliable blocks
Once all the support information for each platform is generated, low and high coverage regions are merged and good supporting regions of each platform are obtained.

This is automatically generated from the above `*/ast.sh` scripts.
 
```
# Get too-low and too-high coverage low support regions by merging blocks less than 100bp away
bedtools subtract -a asm.bed -b ${platform}.bed | bedtools merge -d 100 -i - > $platform.low_high.bed
 
# Trim off regions overlapping 1kb
bedtools subtract -a $platform.low_high.bed -b asm.ends.bed -A > $platform.low_high.trim1k.bed
 
# Get supporting region, using the low_high.bed to exclude <100bp blocks in between other blocks
bedtools subtract -a asm.bed -b $platform.low_high.bed > $platform.support.bed
```

Now, we accumulate the supporting evidence of all platforms and obtain supporting regions where >= 2 platforms agree.

```
Usage: reliable_block.sh <genome> <genomesize.map> [chromosome_assignments.csv]

Run after asset.sh has fully finished its submitted jobs.
Required files: asm.bed, asm.ends.bed, gaps.bed, and */*.support.bed

   <genome>: top level folder of asset
   <genomesize.map>: genome <tab> genome_size (bp)
   [chromosome_assignments.csv]: Scaffold <tab> Assigned_Chromosome <tab> Localized(y/n)
```

`reliable_blocks.sh` is running:
``` 
# Accumulate supports
$asset/bin/acc gaps.bed */*.support.bed > acc.bed 2> acc.log
 
# Merge to get reliable blocks
awk '$4>1' acc.bed | bedtools merge -i - > acc.gt2.mrg.bed
 
# Get low support regions by merging blocks <100bp apart
bedtools subtract -a asm.bed -b acc.gt2.mrg.bed | bedtools merge -d 100 -i - > low_support.bed
 
# Get the final support region as reliable blocks
bedtools subtract -a asm.bed -b low_support.bed > reliable.bed
 
# Exclude low supports in <1kb scaffold boundaries for excluding end-scaffold effects
bedtools subtract -A -a low_support.bed -b asm.ends.bed > low_support.trim1k.bed
```

When `genomesize.map` and `chromosome_assignments.csv` is provided, this script will generate reliable block NG50 stats.
Make sure the top level folder of asset and the genome name are identical in `genomesize.map`.

Example `genomesize.map` used for the first 16 VGP genomes:
```
aRhiBiv1	5067838282
bCalAnn1	1116472572
bStrHab1	1193409022
bTaeGut1	1035611271
bTaeGut2	1009982966
fAnaTes1	662696525
fArcCen1	988048114
fAstCal1	972478608
fCotGob3	608000000
fGouWil2	1182215999
fMasArm1	756753344
mLynCan4	2471525315
mOrnAna1	2128226567
mPhyDis1	2213798723
mRhiFer1	2369916842
rGopEvg1	2684436481
sAmbRad1	3180312983
bTaeGut2.mat	1009982966
bTaeGut2.pat	1035611271
```

## Output
* `reliable.bed` : Reliable regions with support by >= 2 platforms
* `low_support.trim1k.bed` : Regions with low support by > 2 platforms
