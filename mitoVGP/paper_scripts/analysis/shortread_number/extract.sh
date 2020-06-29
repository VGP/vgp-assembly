#!/bin/bash

mkdir -p ${2}

cd ${2}

aws s3 cp --recursive --include="*" --exclude="*_s.fq" s3://genomeark/species/${1}/${2}/assembly_MT_rockefeller/intermediates/bowtie2_round1/fq/ .

fq1=$(wc -l aligned_${2}_all_1.fq | awk '{print $1}')
fq2=$(wc -l aligned_${2}_all_2.fq | awk '{print $1}')

printf "%s\t%s\t%s\n" ${2} $(( fq1 / 4 )) $(( fq2 / 4 ))  > ${2}.cov
