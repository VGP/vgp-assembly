#!/bin/bash

while read -r VGP_ID GEN_ID; do

sh ./repeats.sh $VGP_ID $GEN_ID

done < list.txt
