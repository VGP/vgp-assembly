## Parsing blast tabular output. 

## test file with two scaffolds (one a true mito contaminant) - scaff_52_741_mito_blast_bMelUnd1.tsv 

import csv
import sys
import pandas as pd
from sys import argv

csv.field_size_limit(sys.maxsize)

## Open tabular format output from blast of scaffolds against NCBI mitochondrial db
tabfile = []
with open(argv[1],'r') as file:
    file = csv.reader(file, delimiter = '\t')
    for line in file: 
        tabfile.append(line)

## put in some sort of ... if sum alignment lengths < total scaff length, do not bother checking for coverage 
dffile = pd.DataFrame(tabfile, columns = ['qseqid','sseqid','qlen','length','qcovhsp','eval','qstart','qend','qcovs'])

## list of all the unique scaffolds named in the blast output 
uniqScaffs = dffile.qseqid.unique()

highCovReportLines = []
highCovScaffs = set()

for uniqScaff in uniqScaffs:

    rows = dffile.loc[dffile['qseqid'] == uniqScaff] ## isolate rows with a particular scaffold name 

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
        
        if percentCov > 98: 
            highCovReport = [uniqScaff,uniqAcc,totalScaffLength,coverage,percentCov]
            highCovReportLines.append(highCovReport)
            highCovScaffs.add(uniqScaff)
        # else:
        #     print ("No scaffolds were identified as mitochondrial contaminants")

        
highCovReports = (pd.DataFrame(highCovReportLines, columns=['Accession_num', 'Scaffold','Scaffold_len','Scaffold_align_cov','Perc_align_cov'])).sort_values('Perc_align_cov', ascending = False)
## just pass a file name for now - can rig something better later. 
highCovReports.to_csv(argv[2], sep="\t")

[print (scaff) for scaff in highCovScaffs] 





        


