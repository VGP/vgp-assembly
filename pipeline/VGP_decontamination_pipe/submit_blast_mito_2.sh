#!/bin/bash 

##Provide an input file name (the query sequences) to the first position and the desired name of the outputfile to the second 
##This is the same as the original submit_blast_mito.sh except the num_alignments and coverage parameters have been removed so I can
##get a look at the whole output. 

inputfilename=$1
outfilename=$2.tsv
outreportname=$3
scaffList=$4

blastn -query $inputfilename -db mito_blast_db -num_threads 32 \
-outfmt "6 qseqid sseqid qlen length qcovhsp evalue qstart qend qcovs" -out $outfilename

echo "python3 parse_mito_blast.py $outfilename $outreportname | tr A-Z a-z | tee -a $scaffList"
python3 parse_mito_blast.py $outfilename $outreportname | tr A-Z a-z | tee -a $scaffList

# echo "grep -i 'scaffold' $outreportname | tr A-Z a-z | awk '{print \$2}' | awk 'NR!=1 {print}'| tee -a $scaffList"
# grep -i 'scaffold' $outreportname | tr A-Z a-z | awk '{print $2}' | awk 'NR!=1 {print}'| tee -a $scaffList




