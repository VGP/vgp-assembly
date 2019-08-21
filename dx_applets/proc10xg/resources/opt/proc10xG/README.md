# proc10xG set of python scripts

A series of python scripts to process data generated using the 10x genomics DNA system. scripts are designed to extract and trim the reads of gem barcode information and primer sequence respectively. Compare the gem barcodes to a whitelist (allowing for 1 mismatch) and annotating the read as having a barcode which MATCHES a barcode int the whitelist, contains MISMATCH1, is AMBIGUOUS (with 1 edit distance, matches multiple whitelisted barcodes), or UNKNOWN does not match any whitelisted barcodes at edit distance 1.

Scripts ready for use
* process_10xReads.py - process fastq files generated from bcl2fastq, longranger mkfastq, or supernova mkfastq
* samConcat2Tag.py - extract the FASTA/FASTQ comment appended to SAM output from bwa mem -C and generates 10x genomics sam tags
* filter_10xReads.py - Filters 10x fastq file (processed with process_10xReads.py) by barcode status and/or barcode reads depth or barcode list. (Plan to also support sam/bam input/output)
* regen_10xReads.py - Returns reads fastq file (processed with process_10xReads.py) to 'original' for suitable for input into longranger or supernova (eg after filtering). (Plan to also support sam input)

Scripts in progress, not ready for use
* profile_mapping.py - profile the gem barcode alignments
* process_mapping.py - map + remap ambiguous alignment using gem barcode to identify correct placement

## proc10xGenomics.py, process raw 10x genomic reads (fastq files)

process fastq files generated from bcl2fastq, longranger mkfastq, or supernova mkfastq.
1. extract gem barcode (default: first 16bp of read one), from both sequence and quality
1. trim primer from read (default: next 7bp of read one), from both sequence and quality
1. compare extracted barcode sequence to the whitelist of barcodes, whitelist is expected
	to be in directory barcodes, relative the python script.
1. label read status as
	1. MATCH - perfect match to a whitelist barcode
	1. MISMATCH1 - edit distance of 1 away from a whitelist barcode
	1. AMBIGUOUS - edit distance of 1 away from multiple whitelisted barcodes
	1. UNKNOWN - greater than an edit distance of 1 away from a whitelist barcode
1. annotate reads by appending status, barcode and trimmed sequence to read ID and output

### Usage
	usage: process_10xReads.py [-h] [--version] [-o OUTPUT_DIR] [-a] [-i]
	                           [-b BCTRIM] [-t TRIM] [-g] [--quiet]
	                           [-1 read1 [read1 ...]] [-2 read2 [read2 ...]]

	process_10xReads.py, to process raw fastq files extracting gem barcodes and
	comparing to a white list

	optional arguments:
	  -h, --help            show this help message and exit
	  --version             show program's version number and exit
	  -o OUTPUT_DIR, --output OUTPUT_DIR
	                        Directory + prefix to output reads, [default: stdout]
	  -a, --all             output all reads, not just those with valid gem
	                        barcode, STATUS will be UNKNOWN, or AMBIGUOUS
	                        [default: False]
	  -i                    output in interleaved format, if -o stdout,
	                        interleaved will be chosen automatically [default:
	                        False]
	  -b BCTRIM, --bctrim BCTRIM
	                        trim gem barcode [default: 16]
	  -t TRIM, --trim TRIM  trim additional bases after the gem barcode [default: 7]
	  -g, --nogzip          do not gzip the output, ignored if output is stdout
	  --quiet               turn off verbose output

	Inputs:
	  10x fastq files to input

	  -1 read1 [read1 ...], --read1 read1 [read1 ...]
	                        read1 of a pair, multiple files can be specified
	                        separated by comma
	  -2 read2 [read2 ...], --read2 read2 [read2 ...]
	                        read2 of a pair, multiple files can be specified
	                        separated by comma

	For questions or comments, please contact Matt Settles <settles@ucdavis.edu>
	process_10xReads.py version: 0.0.2

### Output
Reads are output, where the read ID line is annotated with what was extracted in the form:

WHITELIST.GEM.BC:ORIGINAL.READID 1:N:0:LIBRARY.BC:STATUS_GEM.BC_GEM.BC.QUAL_TRIM.SEQ_TRIM.SEQ.QL

Note: 10x barcode status code, gem sequence/qual and trimmed seq/qual are delimited with a _ rather than : because quality values contain : character.

reads can be output as fastq read1 and fastq read 2 in standard format file, or in interleaved format where read 2 follows read 1 in a single file, this facilitates streaming.

#### additional output to standard error, when verbose is on
When verbose option is turned on (default),

Upon opening a file  
PROCESS FILES	reads1_filename	read2_filenames

Upon completion of writing to a file  
PROCESS FILES	Wrote N reads to output

after every 250,000 reads and the final read the following is printed to stdout  
PROCESS READS	reads analyzed:X|reads/sec:X|barcodes:X|median_reads/barcode:X

detailing the applications progress

and at the end of processing  
PROCESS BARCODE	MATCH: X (X%)  
PROCESS BARCODE	MISMATCH1: X (X%)  
PROCESS BARCODE	AMBIGUOUS: X (X%)  
PROCESS BARCODE	UNKNOWN: X (X%)  

These lines can be grepped out of a stdout file, or straight from the output stream

#### whitelisted barcode count

A whitelisted barcode counts file is produced ([output]_barcodes.txt) containing two columns, the barcode sequence and the number of reads assigned to that barcode. Only barcodes found in the whitelist are output

example:  
TGTACGAGTCGGCTAC	3  
CAACCAAGTTACCGAT	1  
CGAAGCCAGAGGGAAT	1  
TCACGCTCACACTCGG	2  
TCGCGTTTCCAGTACA	3  
TATCTACAGTCGTTTG	1  
CTTAATCAGCCATAAA	1  
TTGCCGTGTTAGTGGG	2  

## samConcat2Tag.py, turn gem barcode and trimmed seq to tags

Process a sam formatted file generated from bwa mem after preprocessing reads with proc10xgenomics.py.
bwa mem with the -C option appends the read comment to the end of the sam mapping line, this is not in
a standard tag format and needs to be further parsed. Most software will error on the sam file with the
appended info without additional processing

SAM tags generated from the script
* ST:Z - Read status
* BX:Z - GEM Barcode ID (with appended sample id '-1'), whitelisted ID
* BC:Z - Library Barcode
* QT:Z - Library Barcode Quality (if Index read not provided then all '!')
* RX:Z - GEM Barcode Sequence
* QX:Z - GEM Barcode Quality
* TR:Z - Primer Sequence
* TQ:Z - Primer Quality

### Usage

	usage: samConcat2Tag.py [-h] [--version] [-o OUTPUT_BASE] [inputsam]

	samConcat2Tag, processes bwa mem sam format where the read comment has been
	appended to the mapping line following process_10xReads.py

	positional arguments:
	  inputsam              Sam file to process [default: stdin]

	optional arguments:
	  -h, --help            show this help message and exit
	  --version             show program's version number and exit
	  -o OUTPUT_BASE, --output_base OUTPUT_BASE
	                        Directory + prefix to output, [default: stdout]

	For questions or comments, please contact Matt Settles <settles@ucdavis.edu>
	samConcat2Tag.py version: 0.0.2

## filter_10xReads.py, filters reads for status and barcodes.

This app analyzes reads already processed by process_10xReads.py for specific barcode status (default: MATCH,MISMATCH1) and/or barcodes. You can filter barcodes in 1 of 2 ways.

1. By Specifying the barcode.txt file output by process_10xReads.py and then a min and/or max read count. e.g keep only barcodes with at least min=100 reads and at most max=400 reads.
2. By specifying a list (single column, 1 valid barcode by line) of barcodes to output.

### Usage

	usage: filter_10xReads.py [-h] [--version] [-s STATUSS) [STATUS(S ...]]
	                          [-m BC_MIN] [-n BC_MAX] [-l] [--stdin]
	                          [-o OUTPUT_DIR] [-i] [-g] [--quiet]
	                          [-B barocode.txt] [-L barocode_list.txt]
	                          [-1 read1 [read1 ...]] [-2 [read2 [read2 ...]]]

	filter_10xReads.py, process read file produced by preprocess_10xReads.py and
	filter for certain STATUS conditions, and/or gem barcodes.

	optional arguments:
	  -h, --help            show this help message and exit
	  --version             show program's version number and exit
	  -s STATUS(S) [STATUS(S) ...], --status STATUS(S) [STATUS(S) ...]
	                        which status condition(s) to filter for, allowable
	                        values are MATCH, MISMATCH1, AMBIGUOUS, and UNKNOWN
	                        [default: ['MATCH', 'MISMATCH1']]
	  -m BC_MIN, --min BC_MIN
	                        Minimum barcode read count to output, -B must be
	                        specified [default: None]
	  -n BC_MAX, --max BC_MAX
	                        Maximum barcode read count to output, -B must be
	                        specified [default: None]
	  -l                    input is in interleaved format [default: False]
	  --stdin               accept input on stdin (must be interleaved)
	  -o OUTPUT_DIR, --output OUTPUT_DIR
	                        Directory + prefix to output reads, [default: stdout]
	  -i                    output in interleaved format, if -o stdout,
	                        interleaved will be chosen automatically [default: False]
	  -g, --nogzip          do not gzip the output, ignored if output is stdout
	  --quiet               turn off verbose output

	Inputs:
	  Preprocessed 10x fastq files, and barcode to input

	  -B barocode.txt, --barcode barcode.txt
	                        barcode.txt file produced by process_10xReads.py, only
	                        needed when using --min, --max flags.
	  -L barocode_list.txt, --list barcode_list.txt
	                        A list of barcodes (single column, 1 barcode per row)
	                        to output.
	  -1 read1 [read1 ...], --read1 read1 [read1 ...]
	                        read1 of a pair (or interleaved format), first
	                        processed by process_10xReads, multiple files can be
	                        specified separated by comma
	  -2 [read2 [read2 ...]], --read2 [read2 [read2 ...]]
	                        read2 of a pair, first processed by process_10xReads,
	                        multiple files can be specified separated by comma

	For questions or comments, please contact Matt Settles <settles@ucdavis.edu>
	filter_10xReads.py version: 0.0.1

## regen_10xReads.py, return reads to origin format

Regen, is intended to undo the extraction of sequence/qual in process_10xReads.py, regenerating fastq files that are suitable for processing in supernova, longranger, etc.

### Usage

	usage: regen_10xReads.py [-h] [--version] [-l] [--stdin] [-o OUTPUT_DIR] [-g]
	                         [--quiet] [-1 [read1 [read1 ...]]]
	                         [-2 [read2 [read2 ...]]]

	process_10xReads.py, to process raw fastq files extracting gem barcodes and
	comparing to a white list

	optional arguments:
	  -h, --help            show this help message and exit
	  --version             show program's version number and exit
	  -l                    input is in interleaved format [default: False]
	  --stdin               accept input on stdin (must be interleaved)
	  -o OUTPUT_DIR, --output OUTPUT_DIR
	                        Directory + prefix to output reads, [default: reads]
	  -g, --nogzip          do not gzip the output, ignored if output is stdout
	  --quiet               turn off verbose output

	Inputs:
	  Preprocessed 10x fastq files

	  -1 [read1 [read1 ...]], --read1 [read1 [read1 ...]]
	                        read1 of a pair (or interleaved format), first
	                        processed by process_10xReads, multiple files can be
	                        specified separated by comma
	  -2 [read2 [read2 ...]], --read2 [read2 [read2 ...]]
	                        read2 of a pair, first processed by process_10xReads,
	                        multiple files can be specified separated by comma

	For questions or comments, please contact Matt Settles <settles@ucdavis.edu>
	regen_10xReads.py version: 0.0.1

## Examples

Process 10x reads outputting to files named testing, output all STATUS codes  
> process_10xReads.py -a -o testing -1 data/CaCon-sm_R1_001.fastq.gz -2 data/CaCon-sm_R2_001.fastq.gz

Filtering reads for status codes MATCH,MISMATCH1 (default) and with barcodes read counts between 100 and 400

> filter_10xReads.py -m 100 -n 400 -B testing_barcodes.txt -o testing_filter -1 testing_R1_001.fastq.gz -2 testing_R2_001.fastq.gz

 After filtering, return reads to original for to reprocess with supernova

 > regen_10xReads.py -o rerun_filtered -1 testing_filtered_R1_001.fastq.gz -2 testing_filtered_R1_001.fastq.gz

 OR, by streaming

> filter_10xReads.py -m 100 -n 400 -B testing_barcodes.txt -1 testing_R1_001.fastq.gz -2 testing_R2_001.fastq.gz | regen_10xReads.py -o rerun_filtered

Map processed reads using bwa mem (-C option to append comment), post process with samConcat2tag.py and
then sort by read id (with gem barcode) and finally output bam file  
> bwa mem -C data/polished_p_ctg.fa testing_R1_001.fastq testing_R2_001.fastq | samConcat2Tag.py | samtools sort -n - | samtools view -h -o mapping.sam -

Map processed reads using bwa mem (-C option to append comment), post process with samConcat2tag.py and
then sort by read id (with gem barcode) and finally output directly to bam file

> bwa mem -t 1 -p -C testdata/polished_p_ctg.fa testing_R1_001.fastq testing_R2_001.fastq | samConcat2Tag.py | samtools sort -n -o mapping.bam -

Process a bwa mem sam file with samConcat2Tag.py to exract comment and create tags then sort by position

> samConcat2Tag.py saved.sam | samtools sort - | samtools view

And all together now,

1. first process reads with process_10xReads.py which extracts the GEM barcode and primer sequence and compares
 the barcode to a white list, marking reads with status. Then appends the status, library barcode, GEM barcode,
 primer sequences and corresponding quality scores to the comment of the read ID and the whitelisted barcode to
 the beginning of the read, in interleaved format
1. Then map to the genome using bwa mem with flags -p (interleaved) and -C (appends comment to the sam file)
1. Next process with samContcat2Tag.py which extracts the appended commend and add tags
1. sort using samtools sort, sorting on reads ID (GEM Barcode)
1. and finally saving output to stderr.out and stdout.out

> process_10xReads.py -a -1 data/CaCon-sm_R1_001.fastq.gz \
  -2 data/CaCon-sm_R2_001.fastq.gz | \
  bwa mem -t 1 -p -C data/polished_p_ctg.fa - | samConcat2Tag.py | samtools sort -n -o mapping.bcmapped.bam - 2> stderr.out > stdout.out
