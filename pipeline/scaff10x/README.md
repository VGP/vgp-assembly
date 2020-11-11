# Scaff10x

## Requirements
* [scaff10x](https://github.com/wtsi-hpag/Scaff10X) installed as `$tools/scaff10x/Scaff10X-4.1/src/scaff10x`

## Workflow
```
Usage: ./_submit_scaff10x_v4.0.sh <genome_id> [dir_path_of_fastq.gz_files] [job_id_to_wait_for]
```

### scaff10x
Scaff10x is run with `scaff10x_v4.1.sh`. It internally runs bwa mem and does the scaffolding.
```
$tools/scaff10x/Scaff10X-4.1/src/scaff10x -nodes $((SLURM_CPUS_PER_TASK-2)) -longread 1 -gap 100 -matrix 2000 -read-s1 12 -read-s2 8 -link 10 -score 20 -edge 50000 -block 50000 -data input.dat asm.fasta $name.scaff10x.fasta
```

### Output
* <genome_id>.scaff10x.fasta
