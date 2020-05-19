# WDL VGP Assembly

This directory contains a WDL implementation of the scaffolding portion 
of VGP Assembly, as well as some of the QC steps.  Falcon assembly and 
Arrow polishing have not been implemented in WDL.

## Quick Start

To run WDL workflows, we recommend using `womtool`
to generate JSON describing workflow inputs, and 
`cromwell` for workflow exectution.  They can be downloaded
[here](https://github.com/broadinstitute/cromwell/releases).
Both require Java to run.  This Quick Start
guide assumes running on a single machine.

```
# generate inputs
java -jar /path/to/womtool-0.50.jar \
    inputs \
    /path/to/vgp-assembly/wdl_pipeline/wdl/scaffold_assembly.wdl \
    >input.json
vimacs input.json # edit file

# run workflow
java -jar /path/to/cromwell-50.jar \
    -i input.json \
     /path/to/vgp-assembly/wdl_pipeline/wdl/scaffold_assembly.wdl
```

### Configuration

To configure the workflow, edit the `input.json` file with absolute paths
to your input data files, and specify the thread count you wish to use
on the server.  For optional fields, the default values should be suitable.

The current release version is `0.1.0`.  It is recommended that you
use this instead of the default value "latest" for the `DOCKER_TAG`
parameter.

### Workflows

The top-level workflows are as follows:

* *scaffold_assembly*:  the simplest scaffolding workflow
and is the recommended point of entry.
* *scaffold_assembly_qc*:  the same scaffolding process, but with BUSCO 
and length stats performed after each step.
* *polish_assembly*: runs [MarginPolish](https://github.com/UCSC-nanopore-cgl/MarginPolish)
on an assembly.  This is not part of the VGP assembly pipeline.
* *shasta_assembly_with_scaffold*: runs [Shasta](https://github.com/chanzuckerberg/shasta)
for assembly, then runs the standard scaffolding workflow.  This is
not part of the VGP assembly pipeline.

Each task is configured with its own workflow.  They can be
configured and run in the same way the top-level workflows are.

## Workflow Design

The current workflow design aims to replicate current functionality and
attempts to minimize any changes to the 
existing project, especially with respect to the scripts in the 
`pipeline/` directory.  To do this, the HPC environment
where the scripts are designed to run needs to be replicated. 
The bulk of this work is
installing specific tool versions and configuring them to be usable by
[CEA-HPC Modules](https://github.com/cea-hpc/modules).  Ideally, 
the operative bash scripts are sourced directly from directories in 
`pipeline/` with no required modification.
Slurm submission scripts (ie `_submit_minimap2.sh`) are rewritten and translated into WDL tasks. 
Each task has a WDL file describing its inputs, outputs, resource 
requirements, and the steps required to organize files and invoke the 
operative scripts.  Each task has a task-specific workflow for testing. 
WDL workflows are the intended entrypoint for users; they import tasks 
individually, enabling modularity and task reuse.


## Repository Structure

All locations refer to the `wdl_pipeline/` directory in the main repository.

### Docker

The code to generate Docker images is found in two directories.
The `docker_base/` directory is the base image which contains shared 
libraries and modules configuration.  The `docker_pipeline/` directories
contain images for specific tasks.  These will generally copy the 
operative scripts into a temp directory, install any task-specific
software or libraries, and then copy the scripts onto the Docker image. 

### WDL

The `wdl/` directory holds all the WDL workflows.  These are separated 
into tasks (in `wdl/tasks`) and workflows (in `wdl/`).  Each task file
generally contains a single application and can be imported into multiple
workflows.  Each task file contains a small workflow that can be used
to run that specific application for single-use or debugging.  However,
generally a user will run an entire workflow such as `scaffold_assembly.wdl`.


## Advanced Setup

### Docker

Currently, the Docker images are hosted on 
[DockerHub](https://hub.docker.com/u/tpesout) under tpesout's account
for ease of development.  To run workflows via Cromwell, the images 
need to be pushed to a remote repository.  If a user prefers to host
their own Docker images for use with the workflow, it is possible to
customize the location used as well as the image tag.  

Each docker directory has a Makefile which can be used to build the 
Docker images.  Running `make` will build the image, and `make push`
will push the image to the above repository.  Invoking the makefile 
while passing in the `repository` argument will replace the current 
repository with the specified value.  This applies to all other variables 
specified at the top of the Makefile: repository, identifier, version, 
and tag.  

```
make repository="quay.io/example" version="1.2.3"
make push repository="quay.io/example" version="1.2.3"
```

Without modification, the Makefiles provided will tag each image with
the specified version and current git commit hash from the 
vgp-assembly repository: "$version--$git_hash".  This can be overridden
by specifying the tag= parameter when running the Makefile.
Additionally, all build images will be tagged with "latest". 
Both tags will be pushed to the remote repository.

The `build_docker.sh` script is provided to build and push all docker
images.  It preserves the functionality where it tags images with the
git hash, but could be modified to produce a tag with only the version.

```
./build_docker.sh [repository] [version]
```

One caveat is that the Bionano Docker image requires a tarball of code
which must already be downloaded on the computer before building the images. 
This must be a specific version of the code, which may be difficult for
some users to acquire for manual compilation.  As such, the bionano 
image image is excluded from this build script, and must be built manually
by entering the image's directory and manually running `make`.

#### Bionano Genomics Statement

Bionano Genomics has agreed to provide the licensed Bionano Solve software 
to enable the VGP consortium to package the VGP pipeline in a container. 
The Bionano Solve software in the VGP pipeline may not be the latest 
version, may not have the most recent security patches and may be 
configured in a way unsupported by Bionano. The VGP fully assumes the 
maintenance and support of this VGP pipeline. Please reach out to 
tpesout@ucsc.edu with issues. For the latest supported Bionano software, 
visit [Bionano Support](https://bionanogenomics.com/support/) .
