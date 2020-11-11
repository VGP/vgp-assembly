# Bionano Hybrid Scaffolding
This script started from Solve 3.2.1 and was later updated to 3.3 (3.3_10252018) for DLE1 (DLS) enzymes.

## Requirements
* [Bionano Solve](https://bionanogenomics.com/support/software-downloads/) installed as `$tools/bionano/Solve3.3_10252018`
* Python 2.7.15
* Perl 5.18.2
* R 3.5.0+

## Workflow
Use the matching `_submit_hybrid_*` script for the enzyme used.
* _submit_hybrid_scaffold_dle1_3.3.sh <name>
* _submit_hybrid_scaffold_dle1.sh <name> (3.2.1)
* _submit_hybrid_scaffold_bspqi_bsssi.sh <name> (3.2.1)
* _submit_hybrid_scaffold_bspqi.sh <name> (3.2.1)

### Inputs
In the same directory, the script assumes we have
* Assembly to scaffold symlinked as `asm.fasta`
* Matching molecule `.cmap` files in one of these forms: `DLE1.cmap` `BSPQI.cmap` `BSSSI.cmap`

### Hybrid Scaffold
Depending on how many enzymes used, one of the following gets submitted.

* hybrid_scaffold_3.3.sh (or hybrid_scaffold.sh for 3.2.1)
* hybrid_scaffold_two.sh

Optimize the `.xml` script accordingly.

### Outputs
* `$name/hybrid_scaffolds/*_NCBI.fasta`: Scaffolded fasta
* `$name/hybrid_scaffolds/*_NOT_SCAFFOLDED.fasta`: Not scaffolded fasta
The “Not scaffolded” file contains the left-over from the scaffolding. Mostly contains short contigs, not scaffolded due to alignment failures.

### Next steps
* Concatenate the two fasta files to one:
```
cat $name/hybrid_scaffolds/*_NCBI.fasta $name/hybrid_scaffolds/*_NOT_SCAFFOLDED.fasta > $asm.fasta
```

* Trim trailing N-bases
This was a bug in Solve 3.2.1 and 3.3.3, where some of the trailing N-bases were kept even in the NCBI.fasta version.
Remove these using this [script](https://github.com/VGP/vgp-assembly/tree/master/pipeline/bionano/trimNs).


### Troubleshooting

* Multithreading
This script is running Solve on a single node, scheduled via slurm.
Solve uses a particular python library for multi threading, which was disabled on our cluster.
To enable multithreading, set `OMP_NUM_THREADS`.

`export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK` is handling this in each `hybrid_scaffold*.sh` script.

* Merge multiple .bnx files
Merge `.bnx` files from the same molecule and generate one `.cmap` file before running this script. This happens in case there were several Bionano runs yielding several `.bnx` files.
```
Merge_bnx.sh <prefix>
    <prefix>: ex. mLynCan4_Saphyr_BspQI for mLynCan4_Saphyr_BspQI_748455.bnx and mLynCan4_Saphyr_BspQI_792415.bnx
    <output>: <prefix>.bnx
```
