#!/bin/bash
# purge_dups 0.0.1

set -x -e -o pipefail
main() {

    echo "Value of ref_fastagz: '$ref_fastagz'"
    echo "Value of raw_reads_pacbio_fastagz: '${raw_reads_pacbio_fastagz[@]}'"
    echo "Value of max_genomesize: '$max_genomesize'"
    echo "Value of core_per_job: '$core_per_job'"
    echo "Value of suffix_primary: '$suffix_primary'"
    echo "Value of suffix_haplotig: '$suffix_haplotig'"

    # The following line(s) use the dx command-line tool to download your file
    # inputs to the local file system using variable names for the filenames. To
    # recover the original filenames, you can use the output of "dx describe
    # "$variable" --name".

    dx cat "$ref_fastagz" | zcat > ref.fa
    mkdir -p pacb_fofn

    busyproc=0
    for i in ${!raw_reads_pacbio_fastagz[@]}; do

        if [[ "${busyproc}" -ge $(($(nproc)/$core_per_job)) ]]; then
            echo Processes hit max
            wait -n 
            busyproc=$((busyproc-1))
        fi
        dx download "${raw_reads_pacbio_fastagz[$i]}" 
        /minimap2-2.17_x64-linux/minimap2 -xmap-pb ~/ref.fa -I "$max_genomesize"G "${raw_reads_pacbio_fastagz_name[$i]}" | gzip -c - > $i.paf.gz  &      
        busyproc=$((busyproc+1))
    done 

    while [[ "${busyproc}" -gt  0 ]]; do
        wait -n # p_id
        busyproc=$((busyproc-1))
    done


    /purge_dups/bin/pbcstat *.paf.gz #(produces PB.base.cov and PB.stat files)
    /purge_dups/bin/calcuts PB.stat > cutoffs 2>calcults.log
    #ls  $PWD/* > ~/pb_input
    /purge_dups/bin/split_fa ref.fa > pri_asm.split
    /minimap2-2.17_x64-linux/minimap2 -xasm5 -DP pri_asm.split pri_asm.split | gzip -c - > pri_asm.split.self.paf.gz


    /purge_dups/bin/purge_dups -2 -T cutoffs -c PB.base.cov pri_asm.split.self.paf.gz > dups.bed 2> purge_dups.log
    /purge_dups/bin/get_seqs dups.bed ref.fa > purged.fa 2> hap.fa 
    basename="$ref_fastagz_prefix"
    basename=${basename//.renamed}
    basename=${basename%_c1}
    basename=${basename%_c2p2}
    mv purged.fa "$basename""$suffix_primary".fasta
    mv hap.fa "$basename""$suffix_haplotig".fasta 
    gzip "$basename""$suffix_primary".fasta 
    gzip "$basename""$suffix_haplotig".fasta 
    #python /purge_dups/scripts/run_purge_dups.py config.txt /purge_dups/src $spid

    primary_fastagz=$(dx upload "$basename""$suffix_primary".fasta.gz --brief)
    dup_fastagz=$(dx upload "$basename""$suffix_haplotig".fasta.gz --brief)


    dx-jobutil-add-output primary_fastagz "$primary_fastagz" --class=file
    dx-jobutil-add-output dup_fastagz "$dup_fastagz" --class=file
    mkdir auxillary_files
    mv purge_dups.log calcults.log *.paf.gz dups.bed PB.base.cov PB.cov.wig pri_asm.split cutoffs pacb_fofn PB.stat auxillary_files
    cd auxillary_files/
    for i in $(ls); do 
        auxillary_file=$(dx upload $i --brief)
        dx-jobutil-add-output auxillary_files "$auxillary_file" --class=array:file
    done
}
