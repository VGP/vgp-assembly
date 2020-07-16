#!/bin/bash

SPECIES=${1}
ID=${2}

rm -fr ${ID}

mkdir -p ${ID}

aws s3 --no-sign-request cp --recursive --exclude="*scrap*" --exclude="*xml*" --exclude="*bai*" --exclude="*pbi*" --include="*.subreads.bam" s3://genomeark/species/${SPECIES}/${ID}/genomic_data/pacbio/ ${ID}

cd ${ID}

ls *bam > bam.ls

awk '{print "samtools view "$1" | awk \x27{sum+=length($10)}END{printf (\"%s\\t%s\\n\", \""$1"\", sum)}\x27 >> base.counts"}' bam.ls | parallel

rm *bam
