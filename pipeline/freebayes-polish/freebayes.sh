#!/bin/bash
set -e

i=$SLURM_ARRAY_TASK_ID
sample=$1

if [ -z $i ]; then
    i=$2
fi

if [[ -z $sample ]]; then
    echo "sample=$sample"
    exit -1
fi

module load freebayes/1.2.0
#### [+] Loading freebayes  1.2.0
module load samtools/1.9
#### [+] Loading samtools 1.9  ... 
# This is co-installed with bcftools

module load bedtools

echo "WARNING: This script is no longer in use. Use freebayes_v1.3 version."

#ref=${sample}_t1
ref=asm
fasta=refdata-$ref/fasta/genome.fa
bam=aligned.bam

LEN=`wc -l $fasta.fai | awk '{print $1}'`

mkdir -p bcf

if [ -e bcf/$i.done ]; then
    echo "freebayes for ${i}th array done."
    exit 0 ## Disable temporarily for testing
fi

# Initialize the bcf.$i.list, which tracks the list of bcf files 'done'.
if [ -e bcf.$i.list ]; then
    rm bcf.$i.list
fi

# Perform freebayes every %100 = $i th line
for j in $(seq $i 100 $LEN )
do
    contig=`sed -n ${j}p $fasta.fai | awk '{print $1}'`
    contig_no_pipe=`echo $contig | sed 's/|/_/g'`
    end=`sed -n ${j}p $fasta.fai | awk '{print $2}'`
    # Skip if bcf is not empty
    if ! [ -s bcf/$contig_no_pipe.bcf ]; then
	echo "\
	freebayes --bam $bam --region=$contig:1-$end --max-coverage 200 -f $fasta | bcftools view --no-version -Ou > bcf/$contig_no_pipe.bcf"
        freebayes --bam $bam --region=$contig:1-$end --max-coverage 200 -f $fasta | bcftools view --no-version -Ou > bcf/$contig_no_pipe.bcf
    else
	echo "Found bcf/$contig_no_pipe.bcf"
    fi
    echo "bcf/$contig_no_pipe.bcf" >> bcf.$i.list
done

out=bcf/$i.bcf

echo "\
bcftools concat -f bcf.$i.list -Ou -o $out"
bcftools concat -f bcf.$i.list -Ou -o $out &&

echo "## Clean up the intermediate bcf files" || exit -1

for bcf_file in $(cat bcf.$i.list)
do
    if [ ! -s $bcf_file ]; then
	echo "Skip rm.."
    else
	echo "\
	rm $bcf_file"
	rm $bcf_file	# Disable this too
    fi
done
echo "\
rm bcf.$i.list"
rm bcf.$i.list

echo "\
touch bcf/$i.done"
touch bcf/$i.done
