#!/usr/bin/env python3

import csv
from posixpath import split
import sys
import pandas as pd
from sys import argv
import argparse
from pickle import TRUE
import re 

csv.field_size_limit(sys.maxsize)

parser = argparse.ArgumentParser(add_help=TRUE)
parser.add_argument("--blastout", help="output from submit_mito_blast2.sh; format 6 csv table generated by blasting the assembly against the mitochondrial database")
args = parser.parse_args()

splitdir = (args.blastout).split("/")
newdir = [s for s in splitdir if "output" in s][0]


# ## Open tabular format output from blast of scaffolds against NCBI mitochondrial db
# tabfile = []
# with open(argv[1],'r') as file:
#     file = csv.reader(file, delimiter = '\t')
#     for line in file: 
#         tabfile.append(line)

# ## put in some sort of ... if sum alignment lengths < total scaff length, do not bother checking for coverage 
# dffile = pd.DataFrame(tabfile, columns = ['qseqid','sseqid','qlen','length','qcovhsp','eval','qstart','qend','qcovs'])

# ## list of all the unique scaffolds named in the blast output 
# uniqScaffs = dffile.qseqid.unique()

# highCovReportLines = []
# highCovScaffs = set()

# for uniqScaff in uniqScaffs:

#     rows = dffile.loc[dffile['qseqid'] == uniqScaff] ## isolate rows with a particular scaffold name 

#     uniqAccs = rows.sseqid.unique()

#     totalScaffLength = (list(rows['qlen']))[0]

#     for uniqAcc in uniqAccs:

#         rowsAcc = rows.loc[rows['sseqid'] == uniqAcc]

#         ## making a list of all the alignment start and end positions 
#         starts = list(rowsAcc['qstart'].astype('int')) 
#         ends = list(rowsAcc['qend'].astype('int'))

#         startsEndsDf = pd.DataFrame(list(zip(starts,ends)),
#                                     columns =['starts','ends'])

#         startsEndsDfSort = startsEndsDf.sort_values('starts')

#         startsSort = list(startsEndsDfSort['starts'])
#         endsSort = list(startsEndsDfSort['ends'])

#         coverage = 0
#         currentpos = 0
#         for i in range(len(startsSort)):
#             alignLength = endsSort[i] - startsSort[i]
#             if startsSort[i] > currentpos:
#                 coverage += alignLength
#                 currentpos = endsSort[i]
#             elif (startsSort[i] < currentpos) and (endsSort[i] > currentpos):
#                 coverage += (endsSort[i] - currentpos)
#                 currentpos = endsSort[i]
        
#         percentCov = (coverage/int(totalScaffLength))*100
        
#         if percentCov > 98: 
#             highCovReport = [uniqScaff,uniqAcc,totalScaffLength,coverage,percentCov]
#             highCovReportLines.append(highCovReport)
#             highCovScaffs.add(uniqScaff)
#         # else:
#         #     print ("No scaffolds were identified as mitochondrial contaminants")

        
# highCovReports = (pd.DataFrame(highCovReportLines, columns=['Accession_num', 'Scaffold','Scaffold_len','Scaffold_align_cov','Perc_align_cov'])).sort_values('Perc_align_cov', ascending = False)
# ## just pass a file name for now - can rig something better later. 
# highCovReports.to_csv(argv[2], sep="\t")

# [print (scaff) for scaff in highCovScaffs] 

## Functions ----

def readfile(infile, dfname):
    # Open tabular format output from blast of scaffolds against NCBI mitochondrial db
    tabfile = []
    with open(infile,'r') as file:
        file = csv.reader(file, delimiter = '\t')
        for line in file: 
            tabfile.append(line)

    ## put in some sort of ... if sum alignment lengths < total scaff length, do not bother checking for coverage 
    dfname = pd.DataFrame(tabfile, columns = ['qseqid', 'sseqid', 'length', 'qstart', 'qend', 'evalue' , 'qlen', 'qcovs','qcovhsp'])
    return (dfname)

def calccov(df):
    highCovReportLines = []
    highCovScaffs = set()

    ## list of all the unique scaffolds named in the blast output 
    uniqScaffs = df.qseqid.unique()

    for uniqScaff in uniqScaffs:

        rows = df.loc[df['qseqid'] == uniqScaff] ## isolate rows with a particular scaffold name 

        uniqAccs = rows.sseqid.unique()

        totalScaffLength = (list(rows['qlen']))[0]

        for uniqAcc in uniqAccs:

            rowsAcc = rows.loc[rows['sseqid'] == uniqAcc]

            ## making a list of all the alignment start and end positions 
            starts = list(rowsAcc['qstart'].astype('int')) 
            ends = list(rowsAcc['qend'].astype('int'))

            startsEndsDf = pd.DataFrame(list(zip(starts,ends)),
                                        columns =['starts','ends'])

            startsEndsDfSort = startsEndsDf.sort_values('starts')

            startsSort = list(startsEndsDfSort['starts'])
            endsSort = list(startsEndsDfSort['ends'])

            coverage = 0
            currentpos = 0
            for i in range(len(startsSort)):
                alignLength = endsSort[i] - startsSort[i]
                if startsSort[i] > currentpos:
                    coverage += alignLength
                    currentpos = endsSort[i]
                elif (startsSort[i] < currentpos) and (endsSort[i] > currentpos):
                    coverage += (endsSort[i] - currentpos)
                    currentpos = endsSort[i]
        
            percentCov = (coverage/int(totalScaffLength))*100
            
            if percentCov > 95: 
                highCovReport = [uniqScaff,uniqAcc,totalScaffLength,coverage,percentCov]
                highCovReportLines.append(highCovReport)
                highCovScaffs.add(uniqScaff + "\n")

    return (highCovReportLines, highCovScaffs)
    
## Calculating coverage ----

blastdf = readfile(args.blastout, dfname = "blastdf")
highCovReportLines = calccov(blastdf)

covreportdir = str(newdir) + "/" + "cov_report.tsv"
mitonamesdir = str(newdir) + "/" + "mito_scaff_names.txt"

## Writing the output to a csv report. 
highCovReports = (pd.DataFrame(highCovReportLines[0], columns=['Scaffold','Acc_number','Scaffold_len','Scaffold_align_cov','Perc_align_cov'])).sort_values('Perc_align_cov', ascending = False)
## just pass a file name for now - can rig something better later. 
highCovReports.to_csv(covreportdir, sep="\t")

with open(mitonamesdir, 'w') as f:
    [f.write(scaff) for scaff in highCovReportLines[1]]
    f.close()





        


