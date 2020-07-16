#!/bin/bash

while read -r SPECIES VGP_ID START END TYPE; do

log=$(echo log_${VGP_ID}_${TYPE}.out)

sbatch --partition=vgl,hpc --cpus-per-task=32 --output=${log} het.sh VGP_dataset/${VGP_ID}* ${VGP_ID} ${SPECIES} ${TYPE} 70 10000000 ${START} ${END}

done < $1
