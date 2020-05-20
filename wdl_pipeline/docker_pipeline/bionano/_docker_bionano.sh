#!/bin/bash

# Fix python/3.6 installation to prepend PATH 
module load python/3.5
python --version
which python

module load python/2.7
python --version
which python

module load perl/5.18.2
perl -v
which perl

module load R
R --version
which R