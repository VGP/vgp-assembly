#!/usr/bin/env python
'''
Copyright 2017 Matt Settles
Created June 8, 2017

bwa mem -C option concatenats the fasta/fastq
CTACATTGTCAAGGGT:E00558:34:HGCJ3ALXX:1:1101:2108:1731   99      000000F 922571  60      127M    =       922961  517     ACTCGGGGAGGTGTTAGCTGCTGCCTCACACATTGGGTTTATAGGCTGAATCTTGTTCTCTTTAGGCTTCCAGAGTTTTCTCAGTTACTATTTCTCCTGTCACATACTCGCTGCTTCTTCTGTCATA JJJJJJ<JJF<7A7FJJJJJJ<JJJAJAJJJFJFFFJ----AJJFJ---7---<FJJ<JF<7FFFJJJFJJAJF-AAFFFFF-AFJF7FF<A--FJJJAF)-7-77<<7--)7)<<--77A7-<--< NM:i:3  MD:Z:74T34A3T13 AS:i:112        XS:i:19 1:N:0:GOOD:CCGATTAA:CTACATTGTCAAGGGT:<AAFFJJFJJFJJJJJ:CCAGTGA:J<FFFJJ

This pulls it out, 9 columns and produces new 10x tags in the bam then writes to out
'''
import sys
import os
import argparse

version_num = "0.0.2"
parser = argparse.ArgumentParser(
    description='samConcat2Tag, processes bwa mem sam format where the read comment has been appended to the mapping line following process_10xReads.py',
    epilog='For questions or comments, please contact Matt Settles <settles@ucdavis.edu>\n%(prog)s version: ' + version_num, add_help=True)
parser.add_argument('--version', action='version', version="%(prog)s version: " + version_num)

# TODO: ADD parameter for sample ID
parser.add_argument('-o', '--output_base', help="Directory + prefix to output, [default: %(default)s]",
                    action="store", type=str, dest="output_base", default="stdout")

parser.add_argument('inputfile', metavar='inputsam', type=str, nargs='?',
                    help='Sam file to process [default: %(default)s]', default="stdin")


args = parser.parse_args()  # uncomment this line for command line support


if args.inputfile == 'stdin':
    # reading from stdin
    insam = sys.stdin
else:
    infile = args.inputfile
    # Start opening input/output files:
    if not os.path.exists(infile):
        sys.exit("Error, can't find input file %s" % infile)
    insam = open(infile, 'r')

base = args.output_base

if base is "stdout":
    out = sys.stdout
else:
    out = open(base + ".sam", 'w')

for line in insam:
    # Comment/header lines start with @
    if line[0] != "@" and len(line.strip().split()) > 2:
        line2 = line.strip().split()
        # Handle TAG:
        # get the final concatenatd tag
        tag = line2[-1]
        if (tag[0:6] in ['1:N:0:', '2:N:0:']):
            tsplit = tag.split(":", 4)
            tsplit2 = tsplit[4].split("_")
            if len(tsplit2) != 5:
                sys.stderr.write("SAMCONCAT\tERROR\tsam file has concatenated info, but its the wrong size")
                sys.exit(1)
            # fixed barcode
            line2 = line2[0:-1]
            line2.extend(["ST:Z:" + tsplit2[0],
                          "BX:Z:" + line2[0].split(":")[0] + "-1",
                          "BC:Z:" + tsplit[3],
                          "QT:Z:" + '!' * len(tsplit[3]),
                          "RX:Z:" + tsplit2[1],
                          "QX:Z:" + tsplit2[2],
                          "TR:Z:" + tsplit2[3],
                          "TQ:Z:" + tsplit2[4]])
            out.write('\t'.join(line2) + '\n')
        else:   # Does not contain a concatenated tag as expected by bwa mem
            out.write('\t'.join(line2) + '\n')
    else:  # Its the header lines, so just put back on the stream/file
        out.write(line)


if base is not None:
    out.close()
