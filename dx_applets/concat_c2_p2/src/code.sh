#!/bin/bash
# concat_two_files 0.0.1


main() {

    echo "Value of input1: '$input1'"
    echo "Value of input2: '$input2'"
    echo "Value of output_extension: '$output_extension'"

    dx download "$input1" -o input1
    dx download "$input2" -o input2

    cat input1 input2 > temp_output.fasta.gz

    basename="$input1_prefix"
    basename=${basename//.renamed}
    basename=${basename%_c2}
    basename=${basename%_p2}
    basename="$basename""$output_extension"

    mv temp_output.fasta.gz "$basename".fasta.gz

    output=$(dx upload "$basename".fasta.gz --brief)

    dx-jobutil-add-output output "$output" --class=file
}
