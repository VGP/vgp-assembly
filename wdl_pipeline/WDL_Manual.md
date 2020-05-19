# WDL VGP Assembly

This directory contains a WDL implementation of the scaffolding portion of VGP Assembly.
Also included are WDL implementations for some of the QC steps, as well as
two other tools: Shasta and MarginPolish.  Falcon assembly and Arrow polishing
have not been implemented in WDL.

## Quick Start

To run WDL workflows, we recommend using `womtool`
to generate JSON describing workflow inputs, and 
`Cromwell` for workflow exectution.  They can be downloaded
[here](https://github.com/broadinstitute/cromwell/releases).
Both require Java to run.  This Quick Start
guide assumes running on a single machine.

```
# generate inputs
java -jar /path/to/womtool-0.50.jar \
    inputs \
    /path/to/vgp-assembly/wdl_pipeline/wdl/workflow.wdl \
    >input.json
vimacs input.json # edit file

# run workflow
java -jar /path/to/cromwell-50.jar \
    -i input.json \
     /path/to/vgp-assembly/wdl_pipeline/wdl/workflow.wdl
```

### Configuration

To configure the workflow, edit the `input.json` file with absolute paths
to your input data files, and specify the thread count you wish to use
on the server.  For optional fields, the default values should be suitable.

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

### Conventions

By convention, task-level variables are specified in camelCase and 
workflow variables are in CAPS_CASE.


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
which must already be downloaded on the computer building the images. 
This must be a specific version of the code, which may be difficult for
some users to acquire for manual compilation.  As such, the bionano 
image must be built manually.

Bionano Genomics has agreed to provide the licensed Bionano Solve software 
to enable the VGP consortium to package the VGP pipeline in a container. 
The Bionano Solve software in the VGP pipeline may not be the latest 
version, may not have the most recent security patches and may be 
configured in a way unsupported by Bionano. The VGP fully assumes the 
maintenance and support of this VGP pipeline. Please reach out to 
tpesout@ucsc.edu with issues. For the latest supported Bionano software, 
visit Bionano Support.

### Workflows

When running `wdltool` to generate inputs, only parameters without 
default values are included.  However, there are sometimes task and 
workflow variables that can be used.  Notably, the DOCKER_REPOSITORY
parameter can be used if Docker images have been pushed to a different
repository, and DOCKER_TAG can be specified to use a specific version
(as opposed to "latest", which may be be from code under development).

Note that when running `cromwell` on a single machine, resource 
specifications in the `runtime` section of a task are ignored.  They
are used when configuring tasks to use a certain number of threads.
