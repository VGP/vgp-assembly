#!/bin/bash

while read id
do

sh trnF.sh ../*Circular*${id}* ../../mitos2/NOVOPlasty/*${id}*/*bed

done<$1
