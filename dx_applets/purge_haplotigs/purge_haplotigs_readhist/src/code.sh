#!/bin/bash
set -e -x -o pipefail

MAX_THREADS=16

main() {

    mv gen_histogram.Rscript /usr/scripts/

    dx download "${bam}" -o "subreads.sorted.bam"
    dx download "${genome}" -o "${genome_name}"

    if [[ "${genome_name}" == *.gz ]]; then

        dx download "${genome}" -o - | zcat | python check_contigs.py > genome.fa

    else
        
        dx download "${genome}" -o - | python check_contigs.py > genome.fa
    fi 

    # limit threads to 16 based on issues testing 
    # with a thread count of 32
    threads=-1
    if [ `nproc` -gt $MAX_THREADS ]; then

        threads=$MAX_THREADS

    else

        threads=`nproc`

    fi

    echo "running purge_haplotigs readhist..."
    purge_haplotigs readhist -b subreads.sorted.bam -g genome.fa -t $threads
    echo "Done!"

    mv subreads.sorted.bam.gencov "${genome_prefix}".gencov

    gencov_file_id=$(dx upload "${genome_prefix}".gencov --brief)
    dx-jobutil-add-output genome_coverage "${gencov_file_id}" --class=file

    coverage_pdf_file_id=$(dx upload "${genome_prefix}".coverage_histo.pdf --brief)
    dx-jobutil-add-output histogram "${coverage_pdf_file_id}" --class=file 

}
