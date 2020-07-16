#!/bin/bash

while read -r VGP_ID GEN_ID; do

VGP=$(seqkit fx2tab  --length --name --header-line --gc ../VGP_dataset/$VGP_ID* | tail -1)
GEN=$(seqkit fx2tab  --length --name --header-line --gc ../Genbank_dataset/$GEN_ID* | tail -1)

printf "%s\t%s\t%s\t%s\t%s\t%s\n" $VGP_ID $GEN_ID $(echo $VGP | awk '{print $2"\t"$3}') $(echo $GEN | awk '{print $2"\t"$3}')

done < list.txt
