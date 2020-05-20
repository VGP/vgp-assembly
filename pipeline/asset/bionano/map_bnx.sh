#!/bin/bash

solve_dir=$tools/bionano/Solve3.3_10252018

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

# Get tech and enzyme
maplist=$2
fl=`sed -n ${i}p $maplist`
fn=`basename $fl`
fn_pref=`echo $fn | cut -d_ -f1`
tech=`echo $fn | cut -d_ -f2`
enzyme=`echo $fn | cut -d_ -f3 | sed 's/.bnx$//g'`
output_dir=bnx_${tech}_${enzyme}
output_file=$output_dir/bionano_${tech}_${enzyme}.bed
if [[ -e "$output_file" ]]; then
	echo "$outpuf_file already exists. Exit."
	exit 0
fi

mkdir -p $output_dir

ln -s $ref
ref=`basename $ref`
ref_prefix=${ref%.*}
ref_cmap=$output_dir/fa2cmap/"$ref_prefix"_"${enzyme^^}"_0kb_0labels.cmap
key_fn=$output_dir/fa2cmap/"$ref_prefix"_"${enzyme^^}"_0kb_0labels_key.txt

if [[ ! -s $ref_cmap ]] || [[ ! -s $key_fn ]]; then
	echo "Run fa2cmap"
	echo "
	perl $solve_dir/HybridScaffold/10252018/scripts/fa2cmap_multi_color.pl -e $enzyme 1 -i $ref -o $output_dir/fa2cmap"
	perl $solve_dir/HybridScaffold/10252018/scripts/fa2cmap_multi_color.pl -e $enzyme 1 -i $ref -o $output_dir/fa2cmap
else
	echo "$ref_cmap and $key_fn found."
fi

cp $fl $output_dir/
query_map=$output_dir/$fn
optn=${tech,,}	# to lowercase

map_path=$output_dir/align/contigs/alignmolvref/merge/exp_refineFinal1
rmap_fn="$map_path"_r.cmap
qmap_fn="$map_path"_q.cmap
xmap_fn="$map_path".xmap

if [[ ! -s $xmap_fn ]]; then
	export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

	echo "
	python $solve_dir/Pipeline/10252018/align_bnx_to_cmap.py --prefix $enzyme --mol $query_map --ref $ref_cmap --ra $solve_dir/RefAligner/7915.7989rel/ --nthreads $SLURM_CPUS_PER_TASK --output $output_dir/align --optArgs $solve_dir/RefAligner/7915.7989rel/optArguments_nonhaplotype_"$optn".xml --pipeline $solve_dir/Pipeline/10252018/"
	python $solve_dir/Pipeline/10252018/align_bnx_to_cmap.py --prefix $enzyme --mol $query_map --ref $ref_cmap --ra $solve_dir/RefAligner/7915.7989rel/ --nthreads $SLURM_CPUS_PER_TASK --output $output_dir/align --optArgs $solve_dir/RefAligner/7915.7989rel/optArguments_nonhaplotype_"$optn".xml --pipeline $solve_dir/Pipeline/10252018/
	echo
fi

echo "
$asset/bin/ast_bion_bnx $rmap_fn $qmap_fn $xmap_fn $key_fn > $output_dir/bionano_"$tech"_"$enzyme".bed 2>ast_bion_bnx_"$tech"_"$enzyme".log
"
$asset/bin/ast_bion_bnx $rmap_fn $qmap_fn $xmap_fn $key_fn > $output_dir/bionano_"$tech"_"$enzyme".bed 2>ast_bion_bnx_"$tech"_"$enzyme".log

