#!/bin/bash
# repeat_mito 0.0.1

set -x -e -o pipefail

main() {

    echo "Value of asm: '$mito_asm'"
    
    dx download "$mito_asm" -o mito_asm

	file_type=$(file mito_asm | grep -Po ': \K[^ ]+' | head -1)
	
	if [[ $file_type == "ASCII" ]]; then
		
		mv mito_asm mito_asm.fasta	
	
	elif [[ $file_type == "gzip" ]]; then
		
		mv  mito_asm mito_asm.gz
		
		gunzip -c mito_asm.gz > mito_asm.fasta
	
	else
		
		echo "Please provide a file in either GZIP or ASCII format."
	
		exit 0
		
	fi
	
	header=$(grep ">" mito_asm.fasta)
	sequence=$(grep -v ">" mito_asm.fasta | tr -d "\n")
	
	echo "Original sequence length: $(echo -n $sequence | wc -c)"
	
	repeated_sequence="$sequence$sequence"

	echo "Repeated sequence length: $(echo -n $repeated_sequence | wc -c)"
	
	echo $header > repeat_mito_asm.fasta
	echo $repeated_sequence >> repeat_mito_asm.fasta
	
	gzip repeat_mito_asm.fasta

	repeated_mito_asm=$(dx upload repeat_mito_asm.fasta.gz --brief)
    dx-jobutil-add-output repeated_mito_asm "$repeated_mito_asm" --class=file


}
