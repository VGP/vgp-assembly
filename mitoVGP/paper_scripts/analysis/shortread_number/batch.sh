#!/bin/bash

while IFS=$'\t' read -r species id
do

sbatch extract.sh ${species} ${id}

done<species.ls
