#!/bin/sh

## A shell pipeline for decontaminating VGP genome assemblies post-scaffolding. 

## Help ########################

Help()
{
    ##command line syntax 
    printf "\n"
    echo "  Syntax for running VGP_decontamination_pipe.sh:"
    printf "\n"
    echo "   sh VGP_decontamination_pipe.sh <scaffolded assembly> <unique ID>"
    echo "   Example: sh VGP_decontamination_pipe.sh bTaeGut2.pri.asm.20211014.fasta bTaeGut2.pri"
    printf "\n"
}

while getopts ":h" option; do 
    case $option in
        h) Help
        exit;;
    esac
done 

## Variables ################

scaffs=$1 ##Input file containing assembled (scaffolded) genome. 
genome_id=$2

name=decontam.$genome_id
newdir=outputs.$genome_id
mkdir -p logs
logs=logs/$genome_id.%A_%a.logs

mkdir $newdir ## Directory for output files
touch $newdir/contam_scaffs_$genome_id.txt ## File to append contaminant scaffold names. 
scaffList=$newdir/contam_scaffs_$genome_id.txt 
touch $scaffList

## Ensure all the nucleotides in the input fasta are in upper case before the file goes through soft-masking with dustmasker. 
echo "Converting any lowercase bases to uppercase."
echo -e "tr a-z A-Z < $scaffs | \
    dustmasker -level 40 -out $newdir/masked_$scaffs -outfmt 'fasta'"
tr a-z A-Z < $scaffs | \
    dustmasker -level 40 -out $newdir/masked_$scaffs -outfmt 'fasta'  

echo "Dustmasking complete - will begin to substitute all soft-masked bases for hardmasking."
tr [:lower:] 'N' < $newdir/masked_$scaffs > $newdir/N_sub_masked_$scaffs
# echo -e "python3 sub_soft_hard_mask_pri_asm.py $newdir/masked_$scaffs $newdir/N_sub_masked_$scaffs "
# ## Substitute soft-masking (lowercase bases) for hard-masking ("N") with in-house py script. 
# python3 sub_soft_hard_mask_pri_asm.py $newdir/masked_$scaffs $newdir/N_sub_masked_$scaffs
tr [:lower:] 'N' < $newdir/masked_$scaffs > $newdir/N_sub_masked_$scaffs


echo -e "Masking steps finished. Passing hard-masked scaffolds to Kraken2 for classification."

name=kraken.$genome_id
cores=32
script=$decontam/RELEASE_THE_KRAKEN3.sh
echo -e "sbatch -p vgl -c $cores --error=$logs --output=$logs $script $newdir $scaffs $scaffList"
sbatch -p vgl -c $cores --error=$logs --output=$logs $script $newdir $scaffs $scaffList | awk '{print $4}' > kraken_jid

name=mitoBlast.$genome_id 
cores=32
script=$decontam/submit_blast_mito_2.sh
echo "Blast assembly against mitochondrial sequence database."
echo -e "sbatch -J $name -p vgl -c $cores --error=$logs --output=$logs \
    $script $newdir/N_sub_masked_$scaffs $newdir/mito_blast_$scaffs $newdir/mito_blast_$scaffs.report $scaffList"
sbatch -p vgl -c $cores --error=$logs --output=$logs \
    $script $newdir/N_sub_masked_$scaffs $newdir/mito_blast_$scaffs $scaffList | awk '{print $4}' > mito_jid

jid1=`cat kraken_jid` ## Creating dependencies to delay removal of scaffolds until kraken and mito-blast are both complete. 
jid2=`cat mito_jid`

script=$decontam/clean_fasta.sh

echo -e "sbatch -p vgl --error=$logs --output=$logs --dependency=afterok:$jid1,$jid2 $script $scaffs $newdir/trimmed_$scaffs $scaffList" 
sbatch -p vgl --error=$logs --output=$logs --dependency=afterok:$jid1,$jid2 \
    $script $scaffs $newdir/trimmed_$scaffs $scaffList $newdir/mito_scaff_names.txt $newdir











