#!/bin/bash
# mitoVGP-re 0.0.1, mitoVGP read extractor
# based on mtDNApipe.sh, mitoVGP pipeline v1.1

set -x -e -o pipefail

main() {

    echo "Value of ID: '$ID'"
    echo "Value of REF: '$REF'"

    ref_name=$REF_name
    ref_prefix=$REF_prefix

    dx download "$REF" -o ${ref_name}

	#define working directory
	W_URL=./out/OUTPUT1/assembly_MT/intermediates

	dw_date=`date "+%Y%m%d-%H%M%S"`;
	log_file="log_mitoVGP-re_$dw_date.txt"

	if ! [[ -e "${W_URL}" ]]; then

		mkdir -p ${W_URL}
		echo "--Project folder created: ${W_URL}"
		echo "--Project folder created: ${W_URL}" > $log_file

	else

		echo "--Folder already present: ${W_URL}. Skipping."
		echo "--Folder already present: ${W_URL}. Skipping." >> $log_file	

	fi

	if ! [[ -e "${W_URL}/log" ]]; then

		mkdir -p ${W_URL}/log
		echo "--Log folder created: ${W_URL}/log"
		echo "--Log folder created: ${W_URL}/log" >> $log_file
		echo "--Analysis started on: $dw_date"
		echo "--Analysis started on: $dw_date" >> $log_file

	else

		echo "--Folder already present: ${W_URL}/log. Skipping."
		echo "--Folder already present: ${W_URL}/log. Skipping." >> $log_file	

	fi

	#copy the user-provided reference mitogenome to the reference folder
	if ! [[ -e "${W_URL}/reference" ]] ; then

		mkdir ${W_URL}/reference
		echo "--Reference folder created: ${W_URL}/reference"
		echo "--Reference folder created: ${W_URL}/reference" >> $log_file

		cp ${ref_name} ${W_URL}/reference/${ref_name}
		echo "--User-provided reference copied into reference folder: ${W_URL}/reference/${ref_name}"
		echo "--User-provided reference copied into reference folder: ${W_URL}/reference/${ref_name}" >> $log_file
	
	else

		echo "--Reference folder already present: ${W_URL}/reference. Skipping."
		echo "--Reference folder already present: ${W_URL}/reference. Skipping." >> $log_file	

		if ! [[ -e "${W_URL}/reference/${ref_name}" ]] ; then

			cp ${ref_name} ${W_URL}/reference/${ref_name}
			echo "--User-provided reference copied into reference folder: ${W_URL}/reference/${ref_name}"
			echo "--User-provided reference copied into reference folder: ${W_URL}/reference/${ref_name}" >> $log_file
	
		else

			echo "--Reference file already present in the reference folder. Skipping."
			echo "--Reference file already present in the reference folder. Skipping." >> $log_file
	
		fi
	
	fi

	if ! [[ -e "${W_URL}/pacbio_bam" ]]; then

		mkdir ${W_URL}/pacbio_bam

	else

		echo "--Folder already present: ${W_URL}/pacbio_bam. Skipping."
		echo "--Folder already present: ${W_URL}/pacbio_bam. Skipping." >> $log_file	

	fi

	#for each Pacbio raw data file do
    for i in ${!READS[@]}
    do

    	bam_name=${READS_name[$i]}
    	bam_prefix=${READS[$i]_prefix}
        dx download "${READS[$i]}" -o ${bam_name}

		if ! [[ ${bam_name} == *scraps* ]] && ! [[ ${bam_name} == *.pbi ]] && [[ ${bam_name} == *.bam ]]; then
			
			/opt/smrtlink/smrtcmds/bin/blasr ${bam_name} ${W_URL}/reference/${ref_name} --bam --out ${W_URL}/pacbio_bam/aligned_${bam_name} --bestn 1 --nproc $(nproc)
			echo "--Reads extracted from file: ${bam_name}. Next."
			echo "--Reads extracted from file: ${bam_name}. Next." >> $log_file
    		
    		dx rm ${bam_name}
 		
 		elif ! [[ -e "${W_URL}/pacbio_bam/aligned_${bam_name}" ]]; then
 		
 			echo "--Reads already aligned: ${bam_name}. Skipping."
 			echo "--Reads already aligned: ${bam_name}. Skipping." >> $log_file 		
 		
 		else
 		
 			echo "--Wrong file format for: ${bam_name}. Skipping."
 			echo "--Wrong file format for: ${bam_name}. Skipping." >> $log_file
 		
		fi
 
    done

	echo "--Read extraction completed."
	echo "--Read extraction completed." >> $log_file

	#organize the files
	#convert to fastq
	if ! [[ -e "${W_URL}/pacbio_MT_extracted_reads" ]]; then

		mkdir ${W_URL}/pacbio_MT_extracted_reads
		echo "--Log folder created: ${W_URL}/pacbio_MT_extracted_reads"
		echo "--Log folder created: ${W_URL}/pacbio_MT_extracted_reads" >> $log_file
		
	for f in ${W_URL}/pacbio_bam/aligned_*.bam; do
	
		filename=$(basename -- "$f")
		filename="${filename%.*}"
		bamToFastq -i $f -fq "${W_URL}/pacbio_MT_extracted_reads/${filename}.fq"
		echo "--Fastq conversion completed for file: $f" 
		echo "--Fastq conversion completed for file: $f" >> $log_file

	done
	
	#merge into a single read file
	if ! [[ -e "${W_URL}/pacbio_MT_extracted_reads/${ID}.fastq" ]]; then
		
	cat ${W_URL}/pacbio_MT_extracted_reads/*.fq > ${W_URL}/pacbio_MT_extracted_reads/${ID}.fastq	
	
	gzip ${W_URL}/pacbio_MT_extracted_reads/${ID}.fastq
	
	echo "--Fastq conversion completed."
	echo "--Fastq conversion completed." >> $log_file	
	echo "--Reads collected into: ${W_URL}/pacbio_MT_extracted_reads/${ID}.fastq.gz" 
	echo "--Reads collected into: ${W_URL}/pacbio_MT_extracted_reads/${ID}.fastq.gz" >> $log_file

	fi

	rm ${W_URL}/pacbio_MT_extracted_reads/*.fq

	fi

	mv $log_file ${W_URL}/log/log_mitoVGP-re_$dw_date.txt
 	
	mkdir -p ./out/OUTPUT2/assembly_MT/intermediates/pacbio_MT_extracted_reads
	mv ${W_URL}/pacbio_MT_extracted_reads/${ID}.fastq.gz ./out/OUTPUT2/assembly_MT/intermediates/pacbio_MT_extracted_reads/${ID}.fastq.gz
	dx-upload-all-outputs --parallel
 	
}
