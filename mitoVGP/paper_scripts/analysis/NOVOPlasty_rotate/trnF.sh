#!/bin/bash

GLN_START=$(grep "trnF" $2 | awk '{print $2}')
GLN_END=$(grep "trnF" $2 | awk '{print $2}')
GSIZE=$(awk 'BEGIN {FS="\t"} $0 !~ ">" {sum+=length($0)} END {print sum}' $1)

if (( ${GLN_START} > ${GLN_END} )); then

	printf "\nThe sequence is likely reversed. Generating reverse-complement.\n\n"

	printf "$(sed -n 1p $1)\n$(grep -v ">" $1 | tr -d "\n " | tr "[ATGCUatgcuNnYyRrSsWwKkMmBbDdHhVv]" "[TACGAtacgaNnRrYySsWwMmKkVvHhDdBb]" | rev)" > rv_$(basename -- "${1}")

	tRNAscan-SE rv_${1} -M vert -q > rv_${1}.bed
	
	GLN_START=$(grep "trnF" rv_${1}.bed | awk '{print $3}')
	GLN_END=$(grep "trnF" rv_${1}.bed | awk '{print $4}')
	
	printf "$(sed -n 1p rv_${1})\n$(sed -n 2p rv_${1} | cut -c${GLN_START}-${G_SIZE})$(sed -n 2p rv_${1} | cut -c1-$((${GLN_START} - 1)))" > rotated_${1}

else

	printf "$(sed -n 1p ${1})\n$(cat $1 | grep -v ">" | tr -d "\n" | cut -c${GLN_START}-${G_SIZE})$(cat $1 | grep -v ">" | tr -d "\n " | cut -c1-$((${GLN_START} - 1)))" > rotated_$(basename -- "${1}")

fi

printf "\n" >> rotated_$(basename -- "${1}")
