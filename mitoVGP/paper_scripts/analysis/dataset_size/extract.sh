#!/bin/bash

while IFS=$'\t' read -r species id
do

awk -v ID=${id} '{sum+=$2}END{print ID"\t"sum}' ${id}/base.counts

done<$1
