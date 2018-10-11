#!/bin/bash

if [ -z $1 ]; then
	echo "Usage: ./plot.sh <genomeid>"
	exit -1
fi


module load mash
mash paste combined.msh *.msh
mash dist -t combined.msh combined.msh > combined.tbl

head -n 1 combined.tbl |awk '{for (i=2; i <=NF; i++) print $i}' |awk -F "/" '{print $NF}' |sed s/.subreads.fast[aq].gz//g |sed s/.fast[aq].gz//g |sed s/.fast[aq]//g > key

module load R
Rscript $VGP_PIPELINE/mash/plot.R $1
