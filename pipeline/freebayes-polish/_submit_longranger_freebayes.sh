#! /bin/bash

if [[ -z $1 ]] ; then
        echo "Usage: ./_submit_longranger_freebayes.sh <genome> [ref.fasta]"
        echo "<ref.fasta> will be linked as asm.fasta."
        echo "Assumes we have 10x reads located under /data/rhiea/<genome>/genomic_data/10x/."
        exit -1
fi

if ! [ -e asm.fasta ]; then
        ln -s $2 asm.fasta
fi

sample=$1

mkdir -p logs

cpus=2
mem=12g
name=$sample.longrgr
script=$VGP_PIPELINE/longranger/longranger.sh
args=$sample
walltime=4-0
log=$PWD/logs/$name.%A_%a.log

# launch longranger align
if ! [ -e $sample/outs/possorted_bam.bam ]; then
	echo "\
	sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=norm --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > longranger_jid
	wait_for="--dependency=afterok:`cat longranger_jid`"
fi

if ! [ -e aligned.bam 2> /dev/null ]; then	# symlink regardless the actual file exists or not
	ln -s $sample/outs/possorted_bam.bam aligned.bam 2> /dev/null
	ln -s $sample/outs/possorted_bam.bam.bai aligned.bam.bai 2> /dev/null
	ln -s $sample/outs/summary.csv	2> /dev/null
fi

cpus=4
mem=12g
name=$1.freebayes
script=$VGP_PIPELINE/freebayes-polish/freebayes_v1.3.sh
args=$sample
walltime=3-0
log=logs/$name.%A_%a.log

mkdir -p bcf

fb_done=`ls bcf/*.done | wc -l | awk '{print $1}'`
if ! [[ $fb_done -eq 100 ]]; then
	# Submit job arrays for every 100th contig
	echo "\
	sbatch --partition=norm --array=1-100 -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
	sbatch --partition=norm --array=1-100 -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > freebayes_jid
	wait_for="--dependency=afterok:`cat freebayes_jid`"
fi

cpus=2
mem=4g
name=$sample.consensus
script=$VGP_PIPELINE/freebayes-polish/consensus.sh
args=$sample
walltime=1-0
log=logs/$name.%A_%a.log

echo "\
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > consensus_jid
wait_for="--dependency=afterok:`cat consensus_jid`"

cpus=4
mem=4g
name=$sample.genomecov
script=$VGP_PIPELINE/qv/genomecov.sh
args=$sample
walltime=1-0
log=logs/$name.%A_%a.log

if ! [ -z $3 ]; then
        wait_for="--dependency=afterok:$3"
fi
echo "\
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args"
sbatch --partition=norm -D $PWD $wait_for --cpus-per-task=$cpus --job-name=$name --mem=$mem --time=$walltime --error=$log --output=$log $script $args > genomecov_jid
