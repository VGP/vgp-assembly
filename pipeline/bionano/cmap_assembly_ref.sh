#!/bin/bash

module load python
#### Python 2.7.15 :: Anaconda custom (64-bit)
module load perl/5.18.2
#### Loading Perl 5.18.2  ... 
module load R
#### Loading gcc  7.2.0  ... 
#### Loading GSL 2.4 for GCC 7.2.0 ... 
#### Loading openmpi 3.0.0  for GCC 7.2.0 
#### Loading R 3.5.0_build2

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

BNX=$1
OUT=$2
bionano_path=$tools/bionano/Solve3.2.1_04122018
RUN=$bionano_path/Pipeline/04122018/pipelineCL.py
REFALIGNER=$bionano_path/RefAligner/7437.7523rel/avx/
# Choose the right xml setting
ARG=$3	#$bionano_path/RefAligner/7437.7523rel/optArguments_nonhaplotype_saphyr.xml

python ${RUN} --help
echo

echo "\
python ${RUN} -T $SLURM_CPUS_PER_TASK -i 5 -b ${BNX} -l ${OUT} -t ${REFALIGNER} -a ${ARG} -V 0 -A -m -r ref.cmap"
python ${RUN} -T $SLURM_CPUS_PER_TASK -i 5 -b ${BNX} -l ${OUT} -t ${REFALIGNER} -a ${ARG} -V 0 -A -m -r ref.cmap

