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
#### [+] Loading samtools 1.8  ... 
# This is co-installed with bcftools

module load bedtools


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
	if [ ! -s $contig.bam ]; then
		echo "=== Extract $contig.bam ==="
	        echo "\
		samtools view -hb -@$SLURM_CPUS_PER_TASK -F0x4 $bam $contig:1-$end > $contig.bam"
		samtools view -hb -@$SLURM_CPUS_PER_TASK -F0x4 $bam $contig:1-$end > $contig.bam
		echo
	fi

	if [ ! -s $contig.bed ]; then
		echo "=== Extract $contig.bed: only regions with <50x, larger than 120bp will be saved ==="
		## genomecov is extremely slow even on extracted bam. all tigs in the bam header are reported. 
		## 0-cov region are included (-d option) to avoid too much splitting.
		echo "\
		bedtools genomecov -d -ibam $contig.bam | awk -v contig=$contig '\$1==contig' | java -jar -Xmx1g $VGP_PIPELINE/freebayes-polish/covToRegionBed.jar 150 - | awk '\$3-\$2 > 100 {print \$1\"\t\"\$2+10\"\t\"\$3-10}' - > $contig.bed"
		bedtools genomecov -d -ibam $contig.bam | awk -v contig=$contig '$1==contig' | java -jar -Xmx1g $VGP_PIPELINE/freebayes-polish/covToRegionBed.jar 150 - | awk '$3-$2 > 100 {print $1"\t"$2+10"\t"$3-10}' - > $contig.bed
	fi
        
	if [ -s $contig.bed ]; then
		echo "=== Finally, run freebayes ==="
		echo "\
	        samtools view -F0x4 -hb -L $contig.bed $bam | freebayes -c -f $fasta | bcftools view --no-version -Ou > bcf/$contig_no_pipe.bcf"
        	samtools view -F0x4 -hb -L $contig.bed $bam | freebayes -c -f $fasta | bcftools view --no-version -Ou > bcf/$contig_no_pipe.bcf
		echo

		echo "bcf/$contig_no_pipe.bcf ready!"
		echo "bcf/$contig_no_pipe.bcf" >> bcf.$i.list
		echo
	else
		echo "No region to polish in $contig. Too repetitive."
		echo
	fi
    else
	echo "Found bcf/$contig_no_pipe.bcf"
	echo
	echo "bcf/$contig_no_pipe.bcf" >> bcf.$i.list
    fi
    #echo "bcf/$contig_no_pipe.bcf" >> bcf.$i.list
done

out=bcf/$i.bcf

echo "=== Concatenate bcfs from bcf.$i.list ==="
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

for j in $(seq $i 100 $LEN )
do
    contig=`sed -n ${j}p $fasta.fai | awk '{print $1}'`
    rm $contig.bam
done

echo "\
rm bcf.$i.list"
rm bcf.$i.list

echo "\
touch bcf/$i.done"
touch bcf/$i.done
