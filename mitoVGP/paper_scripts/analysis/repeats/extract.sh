#!/bin/bash

for f in $(ls */summary.txt); do printf "%s\t%s\t%s\t%s\t%s\n" $(head -1 $f | awk '{print $1"\t"$2}') $(tail -1 $f); done
