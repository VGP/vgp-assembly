#!/bin/bash
# proc10xG 0.0.1

set -x -e -o pipefail

main() {

    echo "Value of asm: '$input_fastq'"

	dx-download-all-inputs --parallel
	
	mv in/*/*/*.fastq.gz .
	
	for e in ${input_fastq_name[@]}; do
		
		echo $e
		
		if [[ $e == *"R1"* ]]; then
  	
  		echo $e >> fw.ls
  		
  		elif [[ $e == *"R2"* ]]; then
  		
  		echo $e >> rv.ls
  		
  		fi
	
	done
	
	cat fw.ls | sort | uniq > fw_def.ls
	
	cat rv.ls | sort | uniq > rv_def.ls
	
	cat fw_def.ls
	
	cat rv_def.ls
	
	cat fw_def.ls | grep -oe "S[0-9]_[A-Z][0-9]*" > prefix.ls
	
	paste fw_def.ls rv_def.ls prefix.ls | awk '{print "python /opt/proc10xG/process_10xReads.py -a -1 "$1" -2 "$2" -o trimmed_"$3}' #| parallel --gnu -j $(nproc)

	paste fw_def.ls rv_def.ls prefix.ls | awk '{print "python /opt/proc10xG/process_10xReads.py -a -1 "$1" -2 "$2" -o trimmed_"$3}' | parallel --gnu -j $(nproc)
	
	mkdir -p ~/out/output_fastq
	mv trimmed_*.fastq.gz ~/out/output_fastq

	dx-upload-all-outputs --parallel


}
