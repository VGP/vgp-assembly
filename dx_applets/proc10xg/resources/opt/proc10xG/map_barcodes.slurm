#!/bin/bash
#
#SBATCH --job-name=10x # Job name
#SBATCH --nodes=1
#SBATCH --ntasks=20 # Number of cores
#SBATCH --mem=64000 # Memory pool for all cores (see also --mem-per-cpu)
#SBATCH --time=1-00
#SBATCH --partition=gc128 # Partition to submit to
#SBATCH --output=map10x.out # File to which STDOUT will be written
#SBATCH --error=map10x.err # File to which STDERR will be written
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=youremail@ucdavis.edu # Email to which notifications will be sent

start=`date +%s`

hostname

module load bwa/0.7.13
module load samtools/1.4.1

# requires bwa index on fasta and samtools faidx on fasta
rawbase='/share/biocore/projects/Connon_R_UCD/Menidia_10xgenome/Medinia_Male_L/Project_RCBD_L4_Menidia_Male_L'
fastq1=${rawbase}/Menidia_Male_L_S1_L004_R1_001.fastq.gz
fastq2=${rawbase}/Menidia_Male_L_S1_L004_R2_001.fastq.gz
fasta='/share/biocore/projects/Connon_R_UCD/Menidia_10xgenome/Acanthochromis_polyacanthus/GCA_002109545.1_ASM210954v1_genomic.fna'
output='Menidia-10xMapping2Apolyacanthus.bam'

# assumes that proc10xG repository is in the current working directory
call="proc10xG/process_10xReads.py -a \
 -1 ${fastq1} \
 -2 ${fastq2} \
 | bwa mem -t 16 -p -C \
 ${fasta} - \
 | proc10xG/samConcat2Tag.py \
 | samtools sort -m 768M --threads 4 -n \
 | samtools view -hbt ${fasta}.fai \
 -o ${output} -"


echo $call
eval $call

end=`date +%s`

runtime=$((end-start))

echo $runtime
