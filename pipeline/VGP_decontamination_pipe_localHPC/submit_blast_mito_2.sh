#!/bin/bash 

##Provide an input file name (the query sequences) to the first position and the desired name of the outputfile to the second 
##This is the same as the original submit_blast_mito.sh except the num_alignments and coverage parameters have been removed so I can
##get a look at the whole output. 

inputfilename=$1
outfilename=$2.tsv
scaffList=$3

blastn -query $inputfilename -db $decontam/mito_blast_db -num_threads 32 \
-outfmt "6 qseqid sseqid length qstart qend evalue qlen qcovs qcovhsp" -out $outfilename

echo "python3 parse_mito_blast.py $outfilename | tr A-Z a-z | tee -a $scaffList"
python3 $decontam/parse_mito_blast.py --blastout $outfilename 

# echo "grep -i 'scaffold' $outreportname | tr A-Z a-z | awk '{print \$2}' | awk 'NR!=1 {print}'| tee -a $scaffList"
# grep -i 'scaffold' $outreportname | tr A-Z a-z | awk '{print $2}' | awk 'NR!=1 {print}'| tee -a $scaffList




