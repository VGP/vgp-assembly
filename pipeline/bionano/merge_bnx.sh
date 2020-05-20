#!/bin/bash

if [[ -z $1 ]]; then
	echo "Usage: ./merge_bnx.sh <prefix>"
	echo "    <prefix>: ex. mLynCan4_Saphyr_BspQI for mLynCan4_Saphyr_BspQI_748455.bnx and mLynCan4_Saphyr_BspQI_792415.bnx"
	echo "    <output>: <prefix>.bnx"
	exit 0
fi

prefix=$1

ls ${prefix}*.bnx > $prefix.list

module load perl/5.16.3
module load python/2.7
module load R

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
## Uses at most 4 cpus, ~ 48Gb memory

## Merge the two BspQI .bnx files
echo "
$tools/bionano/Solve3.3_10252018/RefAligner/7915.7989rel/RefAligner -if $prefix.list -merge -o $prefix -bnx -stdout -stderr"
$tools/bionano/Solve3.3_10252018/RefAligner/7915.7989rel/RefAligner -if $prefix.list -merge -o $prefix -bnx -stdout -stderr

