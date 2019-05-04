#!/bin/bash
# concat_two_files 0.0.1


main() {

    echo "Value of s3_fastagz: '$s3_fastagz'"
    echo "Value of c2_fastagz: '$c2_fastagz'"
    echo "Value of p2_fastagz: '$p2_fastagz'"
    echo "Value of outputprefix: '$outputprefix'"

    dx download "$s3_fastagz" -o s3_fasta.gz
    dx download "$c2_fastagz" -o c2_fasta.gz
    dx download "$p2_fastagz" -o p2_fasta.gz

    cat c2_fasta.gz p2_fasta.gz > "$outputprefix"_q2.fasta.gz
    cat s3_fasta.gz "$outputprefix"_q2.fasta.gz > "$outputprefix"_s4.fasta.gz

    q2_fastagz=$(dx upload "$outputprefix"_q2.fasta.gz --brief)
    dx-jobutil-add-output q2_fastagz "$q2_fastagz" --class=file
    s4_fastagz=$(dx upload "$outputprefix"_s4.fasta.gz --brief)
    dx-jobutil-add-output s4_fastagz "$s4_fastagz" --class=file
}
