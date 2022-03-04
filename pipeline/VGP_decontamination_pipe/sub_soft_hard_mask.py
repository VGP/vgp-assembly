## Substituting lower case masking for Ns (hardmask).

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sys import argv
import re

lower_infile = argv[1]
N_outfile = argv[2]

masked_fasta = SeqIO.parse(lower_infile,'fasta')

N_sub_masked_seqs = []
for record in masked_fasta:

    ## Substituting lowercase letters for Ns in all sequences.
    seq = str(record.seq)
    newseq = re.sub("[a-z]","N",seq)

    ## Creating new fasta sequence records with the N-substituted strings. 
    newrecord = SeqRecord(
        Seq(newseq),
        id=record.id,
        description=record.description
    )
    
    N_sub_masked_seqs.append(newrecord)

with open(N_outfile, "w") as output_handle:
    SeqIO.write(N_sub_masked_seqs, output_handle, "fasta")




