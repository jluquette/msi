#!/usr/bin/env python

"""
histogram_msitools.py - Make a frequency table of #supporting reads per locus.
Compute some other interesting stats as well.  Stats in brackets would be
interesting but require joining the original repeat database.
  * how many sites have evidence for multiple alleles
  * mapQ histogram
  * strand distribution
  * chrom distribution
  * STR len (in bp) distribution
  * [ STR type distribution: mono, di, tri, tetra ? ]
  * [ STR genomic location distribution: intergenic, intronic, exonic ? ]
"""

import sys
from argparse import ArgumentParser
from collections import defaultdict
from operator import itemgetter
from numpy import array

parser = ArgumentParser()
parser.add_argument('str_summary', metavar='str_summary', type=str,
                    help='STR summary file from msitools')
parser.add_argument('--mapq60', action='store_true',
                    help='Only consider the highest mapping quality reads (mapQ=60)')

args = parser.parse_args()


def print_sorted(title, d):
    print(title + ":")
    for k in sorted(d.iterkeys()):
        print("%s\t%d" % (k, d[k]))

def emptysplit(x, sep):
    """Return [] if x is the empty string, else x.split as usual."""
    return x.split(sep) if x else []


# line format is tab-delimited: chrom, start, end, STR len list, strand list,
# mapq list.
# STR len list is a string of numbers separated by commas.  The numbers
# are the repeat unit lengths in bp; each number represents a single read.
# So a string like 7,7,7,8 means 3 reads support an allele with 7bp of the
# STR and a single read supports an allele with 8bp of the STR.
# The strand and mapq lists are similarly formatted; the strand is either
# + or - and the mapq is a phred scaled mapping quality between 0 and 60.
# Positions in the 3 lists (STR len, strand, mapQ) all correspond to the
# same read.  So if strlen list=7,7,7,8, strand list=+,-,+,- and mapq list
# is 60,20,30,40; then the first supporting read has an STR len of 7bp,
# occurs on the '+' strand and has an alignment quality of 60.
mapq60_only = args.mapq60
n_refstrs = 0
n_suppreads = 0
chrom_hist = defaultdict(int)
replen_hist = defaultdict(int)
nsupp_hist = defaultdict(int)
nalleles_hist = defaultdict(int)
strand_hist = defaultdict(int)
mapq_hist = defaultdict(int)
f = open(args.str_summary, 'r')
f.readline()  # line 1 is the header
for line in f:
    fields = [ x.strip() for x in line.split('\t') ]
    n_refstrs += 1

    # Track alignment quality distribution
    mapqs = [ int(x.strip()) for x in fields[5].split(',') if x ]
    for mq in mapqs:
        mapq_hist[mq] += 1

    # Determine indexes of mapq=60 reads for optional filtering
    valid_reads = [ idx for idx in range(len(mapqs))
                            if not mapq60_only or mapqs[idx] == 60 ]

    # "if x": don't try to int() cast an empty string, just drop it
    strlens = array([ int(x.strip()) for x in fields[3].split(',') if x ])
    strlens = strlens[valid_reads]
    n_suppreads += len(strlens)
    # How many reads support the locus?
    nsupp_hist[len(strlens)] += 1
    # How many different observed alleles are there?
    nalleles_hist[len(set(strlens))] += 1

    # Track locus distribution across chromosomes by number of supporting
    # reads, not by number in the database.
    chrom_hist[fields[0]] += len(strlens)
    # Track the reference size of the repeat locus (NOT the size(s) observed)
    # Again, track by #supp reads, not by presence in database.
    replen_hist[int(fields[2]) - int(fields[1])] += len(strlens)

    # Track strand distribution
    strands = array([ x.strip() for x in fields[4].split(',') if x ])
    strands = strands[valid_reads]
    for s in strands:
        strand_hist[s] += 1

print("Number of reference STR loci: " + str(n_refstrs))
print("Number of supp reads: " + str(n_suppreads))
print_sorted("Distribution of supp reads by chromosome", chrom_hist)
print_sorted("Distribution of supp reads by reference STR len", replen_hist)
print_sorted("Distribution of reference STR loci by supp reads", nsupp_hist)
print_sorted("Distribution of reference STR loci by num of alleles", nalleles_hist)
print_sorted("Strand frequency over supp reads", strand_hist)
print_sorted("Mapping quality distribution over supp reads", mapq_hist)
