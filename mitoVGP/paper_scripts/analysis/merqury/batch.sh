#!/bin/bash

while IFS=$'\t' read -r species id

do

sbatch mito_qv.sh ${species} ${id} VGP_dataset/${id}*.fasta

done<$1
