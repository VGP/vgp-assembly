# Assembly pipeline for the VGP genomes

The pipeline contains wrappers to run tools in each assembly steps.
Each script is written and tested to run on a slurm scheduler.

To begin, add the path to the pipeline in ~/.bash_profile :
```
export VGP_PIPELINE=/path/to/git/vgp-assembly/pipeline
export tools=/path/to/tools/installed
```

Try `./_submit_` in your desired step to see required input files.
Most of the scripts only require a `.fasta` file as input.

The `_submit` are wrappers to submit jobs to launch the actual scripts.



Tools are loaded as modules. Adjust the `module load` part to your own environment.


