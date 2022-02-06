#!/bin/bash

rm -f filenames.ls

for file in ${2}/[^.]*.ktab
do
basename $file >> filenames.ls
done

printf "merging:"
cat filenames.ls

printf "\ 
Fastmerge -T${3} ${2}/${2} $(awk -F. -v path=${2} '{printf path"/"$1" "}' filenames.ls)"
Fastmerge -T${3} ${2}/${2} $(awk -F. -v path=${2} '{printf path"/"$1" "}' filenames.ls)

