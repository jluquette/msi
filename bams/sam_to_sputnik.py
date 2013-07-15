#!/usr/bin/env python

"""
sam_to_sputnik - convert a SAM file to input for Sputnik.

  --remove_mapq0 - This does NOT remove all MAPQ=0 reads--it removes all
                   PAIRS of reads for which both mates have MAPQ=0.
"""

import sys
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('sam', metavar='SAM file', type=str, nargs=1,
                    help='Aligned SAM (use /dev/stdin with samtools view).')
parser.add_argument('--nopair_mapq0', action="store_true",
                    help='Remove read pairs where both mates have MAPQ=0.')
args = parser.parse_args()

# Globally accessed
mate_dict = {}

def print_fasta(fields):
    print(">" + ' '.join(fields[0:9]))
    print(fields[9])


def output_mapq0(fields):
    try:
        mate = mate_dict[fields[0]]
        # Index [4] is the mapping quality
        if mate[4] != '0' and fields[4] != '0':
            print_fasta(fields)
            print_fasta(mate)
        # Always delete the mate: if they were output, great.  If they weren't,
        # then it's a 0,0 MAPQ pair that will never be output.
        del(mate_dict[fields[0]])
    except KeyError:
        mate_dict[fields[0]] = fields


def output_all(fields):
    print_fasta(fields)


print_fn = output_mapq0 if args.nopair_mapq0 else output_all

f = open(args.sam[0], 'r')
for line in f:
    fields = map(str.strip, line.split('\t'))
    print_fn(fields)
