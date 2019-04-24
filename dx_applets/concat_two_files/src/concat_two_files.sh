#!/bin/bash
# concat_two_files 0.0.1


main() {

    echo "Value of input1: '$input1'"
    echo "Value of input2: '$input2'"
    echo "Value of output_name: '$output_name'"

    dx download "$input1" -o input1
    dx download "$input2" -o input2

    cat input1 input2 > "$output_name"

    output=$(dx upload "$output_name" --brief)

    dx-jobutil-add-output output "$output" --class=file
}
