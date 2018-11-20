import sys

def rename_contigs():
    for l in sys.stdin:
        if l.find('.fasta.gz')>0:
            continue
        if l.startswith('>'):
            if l != l.split('|')[0]:
                raise Exception("Input FASTA header contains invalid character '|'")
            l = l.strip().split('|')[0]
        if len(l.strip())>0:
            print l.strip()
  
if __name__ == '__main__':
    rename_contigs()
