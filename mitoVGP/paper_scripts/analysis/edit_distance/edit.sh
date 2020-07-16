#!/bin/bash

mkdir -p ${3}
mkdir -p ${3}/alignments
mkdir -p ${3}/counts

LN=$(tail -1 $1 | tr -d "\n" | wc -c)

JUMP=10

python edit.py ${1} ${2} ${3} ${JUMP}

original=$(( $(wc -l ${4} | awk '{print $1}') / 4 ))
original_gbp=$(awk 'BEGIN{sum=0;}{if(NR%4==2){sum+=length($0);}}END{print sum;}' ${4})

printf "\nstarting reads: %s\n" ${original}
printf "edits\tedit_fraction\treads\tread_fraction\tgbp\tgbp_fraction\n" > summary_${3}.txt
printf "%s\t%s\t%s\t%s\t%s\t%s\n" 0 0 ${original} 100 ${original_gbp} 100 >> summary_${3}.txt

for (( counter=1; $(( ${counter} * ${JUMP})) <= $(( ${2} * ${JUMP})); counter++ ))
do
	
	count=$(( ${counter} * ${JUMP}))

	pbmm2 align --best-n 1 ${3}/*_${count}.fasta ${4} -j ${5} --sort | samtools view -S -b > ${3}/alignments/${count}.bam
	samtools view ${3}/alignments/${count}.bam | awk '{print $1"\t"length($10)}' > ${3}/counts/${count}.ls
	rm ${3}/alignments/${count}.bam
	mapped=$(wc -l ${3}/counts/${count}.ls | awk '{print $1}')
	mapped_gbp=$(awk '{sum+=$2}END{print sum}' ${3}/counts/${count}.ls)
	printf "\nreads in this iteration: %s\n" ${mapped}
	printf "%s\t%s\t%s\t%s\t%s\t%s\n" ${count} $(echo "scale=10; ${count} / ${LN} *100" | bc) ${mapped} $(echo "scale=10; ${mapped} / ${original} * 100" | bc) ${mapped_gbp} $(echo "scale=10; ${mapped_gbp} / ${original_gbp} * 100" | bc) >> summary_${3}.txt

done