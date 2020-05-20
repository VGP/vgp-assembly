#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./busco.sh <assembly.fasta>"
	exit -1
fi

asm=$1
out=`echo $asm | sed 's/.fasta//' | sed 's/.fa//'`
cpus=$SLURM_CPUS_PER_TASK

#source $tools/conda/etc/profile.d/conda.sh
#conda activate busco

module load bamtools
module load blast/2.2.30+	# multithreading ver.
module load hmmer
module load R
module load python/3.6

#1) bamtools/2.5.1   2) blast/2.2.30++   3) hmmer/3.1b2   4) gcc/7.2.0   5) GSL/2.4_gcc-7.2.0   6) openmpi/3.0.0/gcc-7.2.0-pmi2   7) R/3.5   8) python/3.6

export BUSCO_CONFIG_FILE="$tools/BUSCO/busco/config/config.ini"
export PATH="$tools/Augustus/Augustus-3.3.1/bin:$PATH"
export PATH="$tools/Augustus/Augustus-3.3.1/scripts:$PATH"
export AUGUSTUS_CONFIG_PATH="$tools/Augustus/Augustus-3.3.1/config/"

dataset=$tools/BUSCO/dataset

python $tools/BUSCO/busco/scripts/run_BUSCO.py -h
echo

echo "\
python $tools/BUSCO/busco/scripts/run_BUSCO.py -i $asm -o $out -c $cpus -l $dataset/vertebrata_odb9 -m genome -f -t $out.tmp"
python $tools/BUSCO/busco/scripts/run_BUSCO.py -i $asm -o $out -c $cpus -l $dataset/vertebrata_odb9 -m genome -f -t $out.tmp
