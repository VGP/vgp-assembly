#!/bin/bash
set -e -x -o pipefail

check_cutoffs_are_logical() {


    if [ ${high_cutoff} -lt ${low_cutoff} ]; then 

        dx-jobutil-report-error "ERROR: lot cutoff=${low_cutoff} is greater than high cutoff=${high_cutoff}." 

    fi

    if [ ${midpoint} -gt ${high_cutoff} ]; then

        dx-jobutil-report-error "ERROR: midpoint cutoff=${midpoint} is greater than high cutoff=${high_cutoff}."

    fi 


    if [ ${low_cutoff} -gt ${midpoint} ]; then

        dx-jobutil-report-error "ERROR: midpoint cutoff=${midpoint} is less than low cutoff=${low_cutoff}."

    fi


}



main() {

    mv dot_plot.Rscript /usr/scripts

    check_cutoffs_are_logical
   
    dx download "${genome_coverage}" -o "${genome_coverage_name}" 

    cov_stats=""
    if [ "${coverage_stats}" == "" ]; then

        cov_stats="${genome_coverage_prefix}_coverage_stats.csv"

    else
        
        cov_stats="${coverage_stats}"

    fi

    purge_haplotigs contigcov -i "${genome_coverage_name}" -o "${cov_stats}" -l ${low_cutoff} -m ${midpoint} -h ${high_cutoff} -s ${suspect_threshold} -j ${junk_threshold}
    
    dx download "${genome}" -o assembly.fasta.gz
    gunzip assembly.fasta.gz
    mkdir tmp_purge_haplotigs
    mv assembly.fasta tmp_purge_haplotigs
    dx download "${alignment}" -o "${alignment_name}"
    
    if [ "${dotplots}" == "true" ]; then

        purge_haplotigs purge -g tmp_purge_haplotigs/assembly.fasta -c "${cov_stats}" -b "${alignment_name}"  -t `nproc` -a "${aln_cov}" -dotplots

    else

        purge_haplotigs purge -g tmp_purge_haplotigs/assembly.fasta -c "${cov_stats}" -b "${alignment_name}"  -t `nproc` -a "${aln_cov}"

    fi


    gzip "${cov_stats}"
    cov_stats_fid=$(dx upload "${cov_stats}.gz" --brief)
    dx-jobutil-add-output cov_stats "${cov_stats_fid}" --class=file

    gzip curated.fasta    
    curated_assembly_fid=$(dx upload curated.fasta.gz --brief)
    dx-jobutil-add-output curated_assembly "${curated_assembly_fid}" --class=file    

    gzip curated.haplotigs.fasta
    curated_haplotigs_fid=$(dx upload curated.haplotigs.fasta.gz --brief)
    dx-jobutil-add-output curated_haplotigs "${curated_haplotigs_fid}" --class=file    

    if [ -f curated.artefacts.fasta ]; then
        
        gzip curated.artefacts.fasta
        curated_artifacts_fid=$(dx upload curated.artefacts.fasta.gz --brief)
        dx-jobutil-add-output curated_artifacts "${curated_artifacts_fid}" --class=file   
 
    fi

    curated_reassignments_fid=$(dx upload curated.reassignments.tsv --brief)
    dx-jobutil-add-output curated_reassignments "${curated_reassignments_fid}" --class=file
    
    curated_contig_associations_log_fid=$(dx upload curated.contig_associations.log --brief)
    dx-jobutil-add-output curated_contig_associations "${curated_contig_associations_log_fid}" --class=file


    if [ "${dotplots}" == "true" ]; then 
   
        tar -zcf dotplots_reassigned_contigs.tar.gz dotplots_reassigned_contigs
        reassigned_fid=$(dx upload dotplots_reassigned_contigs.tar.gz --brief)
        dx-jobutil-add-output dotplots_reassigned_contigs "${reassigned_fid}" --class=file
        
        tar -zcf dotplots_unassigned_contigs.tar.gz dotplots_unassigned_contigs
        unassigned_fid=$(dx upload dotplots_unassigned_contigs.tar.gz --brief)
        dx-jobutil-add-output dotplots_unassigned_contigs "${unassigned_fid}" --class=file

    fi

}
