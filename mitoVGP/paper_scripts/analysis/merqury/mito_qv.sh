#!/bin/bash

mkdir -p ${2}

cd ${2}

aws s3 cp --recursive --include="*.fastq.gz" --exclude="*I1*" s3://genomeark/species/${1}/${2}/genomic_data/10x/ genomic_data/10x/

cd genomic_data/10x/

ls *R1*fastq.gz > R1.fofn
ls *R2*fastq.gz > R2.fofn

sbatch /merqury/_submit_build_10x.sh 31 R1.fofn R2.fofn ${2} mem=F

wait_file() {
  local file="$1"; shift

  until [ -f $file ] ; do sleep 300; done
  
}

wait_file ${2}.k31.hist

rm *fastq.gz

meryl greater-than 100 ${2}.k31.meryl output ${2}.k31.gt100.meryl

rm -r ${2}.k31.meryl

cd ../..

ln -s genomic_data/10x/${2}.k31.gt100.meryl
ln -s ../${3}
ass=$(basename ${3})

/merqury/_submit_merqury.sh ${2}.k31.gt100.meryl $ass ${2}