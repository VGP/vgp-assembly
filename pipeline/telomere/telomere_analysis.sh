#!/bin/bash

genome=$1

if [[ -z $genome ]]; then
    echo "Usage: ./telomere.sh <genome_id> <threshold> <ends> <asm.fa> [asm.csv]"
    echo "    <threshold>: threshold to define telomere window. Define from 0.0 to 1.0. Recommended: 0.4"
    echo "    <ends>    : report ends.bed only if a window is found within <ends> bp of scaffold ends"
    echo "    <asm.fa>  : assembly fasta file. .fa or .fasta"
    echo "    <asm.csv> : chromosome assignment, if available. OPTIONAL."
    echo "    <output>  : prefix.windows.<threshold>.bed and"
    echo "                prefix.windows.<threshold>.<ends>.ends.bed"
    exit 0
fi

module load bedtools
module load java/12.0.1
cpus=$SLURM_CPUS_PER_TASK

mkdir -p $genome
cd $genome

if [[ -z $cpus ]]; then
    cpus=1
fi

threshold=$2
ends=$3
asm=$4

echo "genome: $genome"
echo "threshold: $threshold"
echo "ends: $ends"
echo "asm: $asm"

if [[ ! -z $5 ]]; then
    csv=$5
    echo "csv: $csv"
fi
echo

if [[ -z $asm ]]; then
    echo "Failed to find fasta file."
    exit -1
fi
ln -s ../$asm 2> /dev/null

echo "Received $asm"
prefix=`basename $asm`
prefix=`echo $prefix | sed 's/.fasta$//g' | sed 's/.fa$//g'`

#:<<'END'
echo "
$VGP_PIPELINE/telomere/find_telomere.sh $asm"
$VGP_PIPELINE/telomere/find_telomere.sh $asm
echo ""
#END

echo "
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereWindows $prefix.telomere 99.9 $threshold > $prefix.windows.$threshold"
java -cp $VGP_PIPELINE/telomere/telomere.jar FindTelomereWindows $prefix.telomere 99.9 $threshold > $prefix.windows.$threshold
echo

echo "Merge telomere motifs in 100bp"
cat $prefix.windows.$threshold | awk '{print $2"\t"$(NF-2)"\t"$(NF-1)}' | sed 's/>//g' | bedtools merge -d 100  > $prefix.windows.$threshold.bed
echo
#END

echo "Find those at end of scaffolds, within < $ends"
cat $prefix.lens | awk -v ends=$ends '{if ($2>(ends*2)) {print $1"\t0\t"ends"\n"$1"\t"($NF-ends)"\t"$NF} else {print $1"\t0\t"$NF}}' > asm.ends.bed

ends=`echo $ends | awk '{printf "%.0f", $1/1000}'`"kb"
#:<<'END'
bedtools intersect -wa -a $prefix.windows.$threshold.bed -b asm.ends.bed > $prefix.windows.$threshold.$ends.ends.bed

if [[ -z $csv ]]; then
	echo "No .csv found. Exit."
	exit -1
fi

echo "Intersect assigned chromosomes"
echo "Scaff" > $genome.assigned
awk -F "," '{print $1}' $csv  >> $genome.assigned

#:<<'END'
echo "Telomere signals in chr assigned scaffolds"
java -jar -Xmx1g $VGP_PIPELINE/utils/txtContains.jar $prefix.windows.$threshold.bed $genome.assigned 1 > $prefix.windows.$threshold.chr.bed
java -jar -Xmx1g $VGP_PIPELINE/utils/txtContains.jar $prefix.windows.$threshold.$ends.ends.bed $genome.assigned 1 > $prefix.windows.$threshold.$ends.chr.ends.bed

echo "Num. of telomeres at chromosome ends"
echo -e "Scaff\tStart\tEnd" > $prefix.windows.$threshold.$ends.chr.ends.u.bed
bedtools intersect -u -a asm.ends.bed -b $prefix.windows.$threshold.bed >> $prefix.windows.$threshold.$ends.ends.u.bed
java -jar -Xmx1g $VGP_PIPELINE/utils/txtContains.jar $prefix.windows.$threshold.$ends.ends.u.bed $genome.assigned 1 >> $prefix.windows.$threshold.$ends.chr.ends.u.bed
#END

echo "Attach ends"
echo -e "Scaff\tSize" > $prefix.lens.key
cat $prefix.lens >> $prefix.lens.key
java -jar -Xmx1g $VGP_PIPELINE/utils/txtVlookup.jar $genome.assigned $prefix.lens.key $genome.assigned.size Scaff Scaff Size
#END

echo "Attach BE"
java -jar -Xmx1g $VGP_PIPELINE/utils/txtVlookup.jar $genome.assigned.size $prefix.windows.$threshold.$ends.chr.ends.u.bed $genome.assigned.size.start Scaff Scaff Start
cat $genome.assigned.size.start | awk -v genome=$genome '{if ($NF==0) {print genome"\t"$0"\tB"} else if ($NF=="NA") {print genome"\t"$0"\t0"} else {print genome"\t"$0"\tE"}}' | sort -k3 -nr > $genome.summary

rm $genome.assigned.*

cd ../

