#!/bin/bash

while IFS=$'\t' read -r species id ref

do

sbatch NOVOplasty.sh $species $id refs/*${ref}*

done < $1
