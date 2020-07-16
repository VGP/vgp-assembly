#!/bin/bash

while IFS=$'\t' read -r id
do

sh extract.sh ${id} NOVOplasty/${id}/Contigs_1*.fasta ../../paper/VGP_dataset/${id}*.fasta

done<$1
