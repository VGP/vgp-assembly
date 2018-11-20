import sys
import argparse


def _parse_args():
    '''Parse the input arguments.'''
    ap = argparse.ArgumentParser(description='Rename contigs')

    ap.add_argument('-p', '--prefix',
                    help='Prefix of contigs.',
                    default='contig',
                    required=False)

    return ap.parse_args()


def rename_contigs(prefix):
    curr_contig = 0
    for l in sys.stdin:
        if l.startswith('>'):
            l = '>{0}_{1:06x}'.format(prefix, curr_contig)
            curr_contig += 1
        print l.strip()

if __name__ == '__main__':
    args = _parse_args()
    rename_contigs(args.prefix)
