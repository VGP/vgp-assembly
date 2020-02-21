#!/bin/bash
# concat_two_files 0.0.1


main() {

    echo "Value of input1: '$input1'"
    echo "Value of input2: '$input2'"
    echo "Value of output_extension: '$output_extension'"

    dx download "$input1" -o input1
    dx download "$input2" -o input2
	
	if ! [[ -z $mito_asm ]]; then

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

		cat input1 input2 repeat_mito_asm.fasta.gz > temp_output.fasta.gz
	
	else

		cat input1 input2 > temp_output.fasta.gz

	fi
	
    basename="$input2_prefix"
    basename=${basename%_q2}
    basename="$basename""$output_extension"

    mv temp_output.fasta.gz "$basename".fasta.gz

    output=$(dx upload "$basename".fasta.gz --brief)

    dx-jobutil-add-output output "$output" --class=file
}
