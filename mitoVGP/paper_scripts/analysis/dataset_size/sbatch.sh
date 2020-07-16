#!/bin/bash

while IFS=$'\t' read -r species id
do

sbatch download.sh ${species} ${id}

done<$1
