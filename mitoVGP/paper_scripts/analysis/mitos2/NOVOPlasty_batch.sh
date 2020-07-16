#!/bin/bash

while read -r ID; do

mkdir -p $3/$ID

sbatch --output $3/$ID/log.out NOVOPlasty_runmitos.sh $1 $3 $ID

done < $2
