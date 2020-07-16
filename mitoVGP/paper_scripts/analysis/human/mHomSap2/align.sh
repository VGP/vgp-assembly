#!/bin/bash

minimap2 -x asm5 -t 24 ${1} ${2} -a --secondary=no | samtools view -S -b -F 4 -F 0x800 | samtools sort > ${3}
