#!/bin/bash
set -e

echo "Usage: ./freebayes_v1.3.1.sh <sample> [array_idx]"
echo -e "\t<sample>: output prefix"
echo -e "\t[array_idx]: 1-100"

i=$SLURM_ARRAY_TASK_ID
sample=$1

if [ -z $i ]; then
    i=$2
fi

if [[ -z $sample ]]; then
    echo "No <sample> given. Exit."
    exit -1
fi

module load freebayes/1.3.1
#### [+] Loading freebayes  1.3.1
module load samtools/1.9
#### [+] Loading samtools 1.8  ... 
# This is co-installed with bcftools

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
    echo $j
    contig=`sed -n ${j}p $fasta.fai | awk '{print $1}'`
    contig_no_pipe=`echo $contig | sed 's/|/_/g'`
    end=`sed -n ${j}p $fasta.fai | awk '{print $2}'`

    # Skip if bcf is not empty
    if ! [ -s bcf/$contig_no_pipe.bcf ]; then
	mean_cov=`tail -n1 summary.csv | awk -F "," '{printf "%.0f\n", $17}'`
	if [[ "$mean_cov" -eq "" ]]; then
		echo "No mean cov found. Set to 50. Regions with >600x (50 x 12) depth of coverage will be ignored."
		mean_cov=50
	fi
	echo "Using $mean_cov for mean coverage"
	echo "\
	freebayes --bam $bam --region=$contig:1-$end --skip-coverage $((mean_cov*12)) -f $fasta | bcftools view --no-version -Ou > bcf/$contig_no_pipe.bcf"
        freebayes --bam $bam --region=$contig:1-$end --skip-coverage $((mean_cov*12)) -f $fasta | bcftools view --no-version -Ou > bcf/$contig_no_pipe.bcf
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
	rm $bcf_file	# Disable this too when debugging
    fi
done
echo "\
rm bcf.$i.list"
rm bcf.$i.list

echo "\
touch bcf/$i.done"
touch bcf/$i.done


