#!/bin/sh

newdir=$1 
scaffs=$2
scaffList=$3

kraken2 --db $decontam/standard_plus_db --conf 0.30 --threads 32 \
	--classified-out $newdir/class_$scaffs --unclassified-out $newdir/unclass_$scaffs --use-names \
	--output $newdir/classification_$scaffs $newdir/N_sub_masked_$scaffs 

echo "grep -i 'scaffold' $newdir/class_$scaffs | sed 's/>//g' | tr A-Z a-z | awk '{print \$1}' | tee -a $scaffList"
grep -i 'scaffold' $newdir/class_$scaffs | sed 's/>//g' | tr A-Z a-z | awk '{print $1}' | tee -a $scaffList
