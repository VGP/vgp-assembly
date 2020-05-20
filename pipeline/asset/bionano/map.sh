#!/bin/bash

solve_dir=/data/Phillippy/tools/bionano/Solve3.2.1_04122018

ref=$1

if [ -z $ref ]; then
	echo "Usage: ./map.sh <ref.fasta> <input.fofn> [line_num]"
	exit -1
fi

if [ -z $2 ]; then
	echo "No input.fofn given. Exit."
	exit -1
fi

i=$SLURM_ARRAY_TASK_ID

if [ -z $i ]; then
	i=$3
fi

if [ -z $i ]; then
	echo "No SLURM_ARRAY_TASK_ID or line_num given. Exit."
	exit -1
fi

module load perl/5.16.3
module load python/2.7
module load R

cmaplist=$2
fl=`sed -n ${i}p $cmaplist`
fn=`basename $fl`
fn_pref=`echo $fn | cut -d_ -f1`
tech=`echo $fn | cut -d_ -f2`
enzyme=`echo $fn | cut -d_ -f3 | sed 's/.cmap$//g'`
output_dir=$enzyme

mkdir -p $output_dir

ln -s $ref
ref=`basename $ref`
ref_prefix=${ref%.*}
ref_cmap=$output_dir/fa2cmap/"$ref_prefix"_"$enzyme"_0Kb_0labels.cmap
key_fn=$output_dir/fa2cmap/"$ref_prefix"_"$enzyme"_0Kb_0labels_key.txt

if [[ ! -s $ref_cmap ]] || [[ ! -s $key_fn ]]; then
	echo "Run fa2cmap"
	echo "
	perl $solve_dir/HybridScaffold/04122018/scripts/fa2cmap.pl -n ${enzyme} -i $ref -o $output_dir"
	perl $solve_dir/HybridScaffold/04122018/scripts/fa2cmap.pl -n ${enzyme} -i $ref -o $output_dir
else
	echo "$ref_cmap and $key_fn found."
fi

cp $fl $output_dir/
query_cmap=$output_dir/$fn
optn=${tech,,}	# to lowercase

echo "
python2 $solve_dir/Pipeline/04122018/runCharacterize.py -t   $solve_dir/RefAligner/7437.7523rel/RefAligner -q $query_cmap -r  $ref_cmap -p $solve_dir/Pipeline/04122018/ -a $solve_dir/RefAligner/7437.7523rel/optArguments_nonhaplotype_"$optn".xml -n $SLURM_CPUS_PER_TASK"
python2 $solve_dir/Pipeline/04122018/runCharacterize.py -t   $solve_dir/RefAligner/7437.7523rel/RefAligner -q $query_cmap -r  $ref_cmap -p $solve_dir/Pipeline/04122018/ -a $solve_dir/RefAligner/7437.7523rel/optArguments_nonhaplotype_"$optn".xml -n $SLURM_CPUS_PER_TASK

echo

map_path=$output_dir/alignref/${fn%.*} 
rmap_fn="$map_path"_r.cmap
qmap_fn="$map_path"_q.cmap
xmap_fn="$map_path".xmap

echo "
$asset/bin/ast_bion $rmap_fn $qmap_fn $xmap_fn $key_fn > $output_dir/bionano_"$tech"_"$enzyme".bed 2>ast_bion_"$tech"_"$enzyme".log
"
$asset/bin/ast_bion $rmap_fn $qmap_fn $xmap_fn $key_fn > $output_dir/bionano_"$tech"_"$enzyme".bed 2>ast_bion_"$tech"_"$enzyme".log
