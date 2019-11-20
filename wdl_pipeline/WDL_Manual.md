# WDL VGP Assembly

An implementation of the VGP Assembly is currently under development.

## Quick Start

To run WDL workflows, we recommend using 
[wdltool](https://github.com/broadinstitute/wdltool)
to generate JSON describing workflow inputs, and 
[Cromwell](https://github.com/broadinstitute/cromwell)
for workflow exectution.  Both require Java to run.  This Quick Start
guide assumes running on a single machine.

```
# link import directory (wdltool needs task imports to be local)
cd /your/work/directory
ln -s /path/to/vgp-assembly/wdl_pipeline/wdl/tasks

# generate inputs
java -jar /path/to/wdltool-0.14.jar \
    inputs \
    /path/to/vgp-assembly/wdl_pipeline/wdl/workflow.wdl \
    >input.json
vimacs input.json # edit file

# run workflow
java -jar /path/to/cromwell-44.jar \
    -i input.json \
     /path/to/vgp-assembly/wdl_pipeline/wdl/workflow.wdl
```

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
Slurm submission scripts are rewritten and translated into WDL tasks. 
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
The `docker_base/` directory is the base image which contains all shared 
libraries and modules configuration.  The `docker_pipeline/` directories
contain images for specific tasks.  These will generally copy the 
operative scripts into a temp directory, install any task-specific
software or libraries, and then copy the scripts onto the Docker image. 

### WDL

The `wdl/` directory holds all the WDL workflows.  These need to import 
individual tasks held in the `wdl/tasks/` directory.  Note that when
running the `cromwell` jar, imported files can be found relative to the
invoked workflow.  However, when running the `wdltool` jar, the imports
need to be relative to the current directory.  The `ln` command in the
Quick Start section handles this.


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

### Workflows

When running `wdltool` to generate inputs, only parameters without 
default values are included.  However, there are sometimes task and 
workflow variables that can be used.  Notably, the DOCKER_REPOSITORY
parameter can be used if Docker images have been pushed to a different
repository, and DOCKER_TAG can be specified to use a specific version
(as opposed to "latest", which may be be from code under development).

Note that when running `cromwell` on a single machine, resource 
specifications in the `runtime` section of a task are ignored.  These
are included in the tasks because they will eventually be used. 
Specifiying THREAD_COUNT options are still useful, as this parameter
is used when configuring tasks to use a certain number of threads.
