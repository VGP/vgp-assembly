#!/bin/bash

set -e

if [ -z $1 ]; then

	echo "use $0 -h for help"
	exit 0
elif [ $1 == "-h" ]; then

	cat << EOF

	Usage: '$0 species id reference genome_size mapper threads'

	mtDNApipe.sh is used to retrieve mitochondrial-like sequences from the raw Pacbio data
	generated in the framework of the Vertebrate Genomes Project and assemble them using Canu.

	It requires the following software (and their dependencies) installed:
	aws-cli/1.16.101, blasr/5.3.2-06c9543 | minimap2/2.17 | pbmm2/1.0.0, bam2fastx, Canu/1.8, blastn/2.7.1+, pbindex

	Sequence retrieval is based on a search by similarity using BLASR alignment.
	Pacbio raw data files are individually downloaded from the Genomeark
	and aligned to a reference genome provided by the user.

	The reference genome can be from the same species if available, or from a
	closely-to-distantly related species.

	The approach is similar to that of Organelle_PBA described in:
	Soorni et al. BMC Genomics (2017) DOI 10.1186/s12864-016-3412-9

	In the second steps reads are used to generate assemblies using Canu assembler,
	usually with default parameters.
	(for rGopEvg1 and mRhiFer1 minOverlapLength=300 correctedErrorRate=0.105 were used, respectively)

	The reference genome provided by the user is then blasted to the contigs generated
	by Canu to identify the putative mitocontig.

	Required arguments are:
	-s the species name (e.g. Calypte_anna)
	-i the VGP species ID (e.g. bCalAnn1)
	-r the reference sequence fasta file
	-g the putative mitogenome size (potentially, that of the reference genome). It does not
	need to be precise. Accepts Canu formatting.
	-t the number of threads
	
	Optional arguments are:
	-l the aligner (blasr|minimap2|pbmm2). Default is pbmm2
	-o the options for Canu

EOF

exit 0

fi

#set options

while getopts ":l:s:i:r:g:i:o:m:t:w:" opt; do

	case $opt in
		l)
			LIST=$OPTARG
			echo "Mode: -l $OPTARG"
			;;
		s)
			SPECIES=$OPTARG
			echo "Species: -s $OPTARG"
			;;
        i)
        	ABBR=$OPTARG
        	echo "Species ID: -i $OPTARG"
            ;;
        r)
			REF=$OPTARG
			echo "Reference: -r $OPTARG"
			;;
		g)
			SIZE=$OPTARG
			echo "Genome size: -g $OPTARG"
            ;;
		c)
            GRID=$OPTARG
			echo "Cluster: -c $OPTARG"
			;;
		o)
            OPTS=$OPTARG
			echo "Canu options: -o $OPTARG"
			;;
		m)
			ALN=$OPTARG
			echo "Aligner: -m $OPTARG" >&2
			;;
		t)
			NPROC=$OPTARG
			echo "Number of threads: -t $OPTARG" >&2
            ;;
		\?)
			echo "ERROR - Invalid option: -$OPTARG" >&2
			exit 1
			;;
	esac

done

if [[  ${GRID} == "SLURM" ]]; then

#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
##SBATCH --ntasks=1
echo Starting at `date`
echo This is job $SLURM_JOB_ID
echo Running on `hostname`

fi

#define working directory
W_URL=${SPECIES}/assembly_MT/intermediates

if ! [[ -e "${W_URL}" ]]; then

	mkdir -p ${W_URL}

fi

#copy the user-provided reference mitogenome to the reference folder
if ! [[ -e "${W_URL}/reference" ]]; then

	mkdir ${W_URL}/reference
	cp ${REF} ${W_URL}/reference/${REF%.*}.fasta

fi

if ! [[ -e "${W_URL}/log" ]]; then

	mkdir ${W_URL}/log

fi

dw_date=`date "+%Y%m%d-%H%M%S"`;

if [[ -z ${LIST} ]]; then

	#record Pacbio raw data files available in the cloud at the time of the analysis
	aws s3 --no-sign-request ls s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/pacbio/ | grep -oP "m.*.subreads.bam" | uniq > ${W_URL}/log/file_list_$dw_date.txt

else

	cat ${LIST} > ${W_URL}/log/file_list_$dw_date.txt

fi

if ! [[ -e "${W_URL}/pacbio_bam" ]]; then

mkdir ${W_URL}/pacbio_bam

fi

if [[ -z  ${ALN} ]]; then
	ALN="pbmm2"
fi

if [[  ${ALN} == "pbmm2" ]] && ! [[ -e ${W_URL}/reference/${REF%.*}.fasta.mmi ]]; then                                                                                                      
	pbmm2 index ${W_URL}/reference/${REF%.*}.fasta ${W_URL}/reference/${REF%.*}.fasta.mmi 

fi

#for each Pacbio raw data file do
while read p; do

if ! [[ $p == *scraps* ]] && ! [[ $p == *.pbi ]] && [[ $p == *.bam ]] && ! [[ -e "${W_URL}/pacbio_bam/aligned_$(basename -- "$p")" ]]; then
	
	if [[ -z ${LIST} ]]; then

		#if vgp mode download
		aws s3 --no-sign-request cp s3://genomeark/species/${SPECIES}/${ABBR}/genomic_data/pacbio/$p ${W_URL}/
		#align

	else
		
		ln -s $p ${W_URL}/$(basename -- "$p")
		p=$(basename -- "$p")
	fi

	if [[  ${ALN} == "minimap2" ]]; then

		pbindex ${W_URL}/${p%.*}.bam

		bam2fastq ${W_URL}/${p%.*}.bam -o ${W_URL}/${p%.*}

		minimap2 -ax map-pb --secondary=no ${W_URL}/reference/${REF%.*}.fasta ${W_URL}/${p%.*}.fastq.gz -t ${NPROC} -R "$(samtools view -H ${W_URL}/${p%.*}.bam | grep "@RG" | tr -d '\n' | sed 's/	/\\t/g')" | samtools view -hbS -F2048 -F4 - > ${W_URL}/pacbio_bam/aligned_${p%.*}.bam
		rm ${W_URL}/${p%.*}.fastq.gz
		rm ${W_URL}/${p%.*}.bam.pbi

	elif [[  ${ALN} == "blasr" ]]; then

		blasr --bestn 1 ${W_URL}/$p ${W_URL}/reference/${REF%.*}.fasta --bam --out ${W_URL}/pacbio_bam/aligned_${p%.*}.bam --nproc ${NPROC}

	elif [[  ${ALN} == "pbmm2" ]]; then
	
		pbmm2 align ${W_URL}/reference/${REF%.*}.fasta.mmi ${W_URL}/$p ${W_URL}/pacbio_bam/aligned_${p%.*}.bam -j ${NPROC} -N 1
		
	else

		echo "mapper unset or unidentified"

	fi
#remove
rm ${W_URL}/${p%.*}.bam

fi

done <${W_URL}/log/file_list_$dw_date.txt

#organize the files

#convert to fastq
if ! [[ -e "${W_URL}/pacbio_MT_extracted_reads" ]]; then

mkdir ${W_URL}/pacbio_MT_extracted_reads

for f in ${W_URL}/pacbio_bam/aligned_*.bam; do
	filename=$(basename -- "$f")
	filename="${filename%.*}"
	if ! [[ -e "${f}.pbi" ]]; then

		pbindex ${f}

	fi

	if ! [[ -e "{W_URL}/pacbio_MT_extracted_reads/${filename}.fastq.gz" ]]; then

		bam2fastq ${f} -o "${W_URL}/pacbio_MT_extracted_reads/${filename}"

	fi
done

#merge into a single read file
if ! [[ -e "${W_URL}/pacbio_MT_extracted_reads/${ABBR}.fastq" ]]; then

zcat ${W_URL}/pacbio_MT_extracted_reads/*.fastq.gz > ${W_URL}/pacbio_MT_extracted_reads/${ABBR}.fastq

fi

rm ${W_URL}/pacbio_MT_extracted_reads/*.fastq.gz

fi

#assemble mtDNA reads with canu
if ! [[ -e "${W_URL}/canu/${ABBR}.contigs.fasta" ]]; then

	CANU=/software/grit/mitoVGP/canu/Linux-amd64/bin/canu

	if ! [[ -z  ${OPTS} ]]; then

		CANU="${CANU} ${OPTS}"

	fi
	
	CANU="${CANU} contigFilter=\"2 0 1.0 0.5 0\" -p ${ABBR} -d ${W_URL}/canu"
	CANU="${CANU} genomeSize=${SIZE}"
#	CANU="${CANU} -gridEngineResourceOption=\"-R\\\"select[mem>${LSF_MEM}] rusage[mem=${LSF_MEM}]\\\" -M${LSF_MEM} -n ${NPROC}\""
	CANU="${CANU} -pacbio-raw ${W_URL}/pacbio_MT_extracted_reads/${ABBR}.fastq"
	
	echo ${CANU}
	eval ${CANU}
fi

wait

if ! [[ -e "${W_URL}/blast" ]] && [[ -e "${W_URL}/canu/${ABBR}.contigs.fasta" ]]; then

mkdir -p ${W_URL}/blast

fi

if [[ -e "${W_URL}/canu/${ABBR}.contigs.fasta" ]]; then

#build blast db
makeblastdb -in ${W_URL}/canu/${ABBR}.contigs.fasta -parse_seqids -dbtype nucl -out ${W_URL}/blast/${ABBR}.db

#search the putative mitocontig using blastn
blastn -query ${W_URL}/reference/${REF%.*}.fasta -db ${W_URL}/blast/${ABBR}.db -out ${W_URL}/blast/${ABBR}.out
blastn -query ${W_URL}/reference/${REF%.*}.fasta -db ${W_URL}/blast/${ABBR}.db -outfmt 6 -out ${W_URL}/blast/${ABBR}.tb
sed -i "1iquery_acc.ver\tsubject_acc.ver\t%_identity\talignment_length\tmismatches\tgap_opens\tq.start\tq.end\ts.start\ts.end\tevalue\tbitscore" ${W_URL}/blast/${ABBR}.tb
cat ${W_URL}/blast/${ABBR}.tb | column -t > ${W_URL}/blast/${ABBR}_results.txt

fi