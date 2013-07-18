#!/usr/bin/env python

"""
metrics_msitools.py - Make a frequency table of #supporting reads per locus.
Compute some other interesting stats as well.  Stats in brackets would be
interesting but require joining the original repeat database.
  * how many sites have evidence for multiple alleles
  * mapQ histogram
  * strand distribution
  * chrom distribution
  * STR len (in bp) distribution
  * difference between observed STR len and ref STR len
    not absolute value: negative values show obs < ref len
  * [ STR unit distribution: mono, di, tri, tetra ? ]
  * [ STR unit multiple distribution: (end-start) / unit ]
  * [ STR genomic location distribution: intergenic, intronic, exonic ? ]

Options for filtering reads:
  * haploid options: --only-x, --only-y, only use X or Y alignments, depending
    on option (both cannot be specified together).  For male samples, these
    should be haploid and thus deviations from expectation should reflect
    experimental errors.
  * maximum mapQ option: --mapq60, only use mapq60 (very high confidence)
    alignments
  * [ minimum repeat unit option: --min-unit (recommended=3, default=1?).  Only
    consider loci where the repeat unit (end-star+1) / unit size is greater
    than the specified value.  This could be useful to remove questionable
    loci, like 2~3 repeat units of tri or tetranucleotide repeats. ]
"""

import sys
from argparse import ArgumentParser
from collections import defaultdict
from operator import itemgetter
from numpy import array

def print_sorted(title, d):
    print(title + ":")
    for k in sorted(d.iterkeys()):
        print("%s\t%d" % (k, d[k]))

def emptysplit(x, sep):
    """Return [] if x is the empty string, else x.split as usual."""
    return x.split(sep) if x else []


parser = ArgumentParser()
parser.add_argument('str_summary', metavar='str_summary', type=str,
                    help='STR summary file from msitools')
parser.add_argument('--mapq60', action='store_true',
                    help='Only consider the highest mapping quality reads (mapQ=60)')
parser.add_argument('--x-only', dest='x_only', action='store_true',
                    help='Only consider reads from the X chromosome')
parser.add_argument('--y-only', dest='y_only', action='store_true',
                    help='Only consider reads from the Y chromosome')
parser.add_argument('--min-units', dest='min_units', metavar='N', default=1,
                    help='Discard reference loci with fewer than N repeat units')

args = parser.parse_args()

if args.x_only and args.y_only:
    print('ERROR: only one of --x-only and --y-only may be specified')
    exit(1)

unit_to_int = { 'mono': 1, 'di': 2, 'tri': 3, 'tetra': 4 }
mapq60_only = args.mapq60
x_only = args.x_only
y_only = args.y_only
min_units = args.min_units
total_refstrs = 0
total_reads = 0
n_refstrs = 0
n_suppreads = 0
loci_x_only_filter = 0
reads_x_only_filter = 0
loci_y_only_filter = 0
reads_y_only_filter = 0
loci_mapq60_filter = 0
reads_mapq60_filter = 0
loci_min_units_filter = 0
reads_min_units_filter = 0
chrom_hist = defaultdict(int)
reflen_hist = defaultdict(int)
reflen_diff_hist = defaultdict(int)
nsupp_hist = defaultdict(int)
nalleles_hist = defaultdict(int)
strand_hist = defaultdict(int)
mapq_hist = defaultdict(int)

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
f = open(args.str_summary, 'r')
f.readline()  # line 1 is the header
for line in f:
    fields = [ x.strip() for x in line.split('\t') ]

    total_refstrs += 1

    mapqs = [ int(x.strip()) for x in fields[5].split(',') if x ]
    total_reads += len(mapqs)

    if x_only and not fields[0].endswith('X'):
        loci_x_only_filter += 1
        reads_x_only_filter += len(mapqs)
        continue

    if y_only and not fields[0].endswith('Y'):
        loci_y_only_filter += 1
        reads_y_only_filter += len(mapqs)
        continue

    reflen = int(fields[2]) - int(fields[1]) + 1  # reference STR length in bp
    # TODO: use this min units filter
    # if float(reflen) / unit_to_int[fields[unit]] < min_units:
    #   loci_min_units_filter += 1
    #   reads_min_units_filter += len(mapqs)
    #   continue


    # Determine indexes of mapq=60 reads for optional filtering
    valid_reads = [ idx for idx in range(len(mapqs))
                            if not mapq60_only or mapqs[idx] == 60 ]

    reads_mapq60_filter += len(mapqs) - len(valid_reads)
    if not valid_reads:
        # If len(mapqs) = 0, there were no reads to begin with, so not an
        # effect of this filter.
        if len(mapqs) > 0:
            loci_mapq60_filter += 1
        continue

    # Track alignment quality distribution
    for mq in mapqs:
        mapq_hist[mq] += 1

    n_refstrs += 1

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
    reflen_hist[reflen] += len(strlens)

    # Track the difference between observed and reference STR lens
    for obslen in strlens:
        reflen_diff_hist[obslen - reflen] += 1

    # Track strand distribution
    strands = array([ x.strip() for x in fields[4].split(',') if x ])
    strands = strands[valid_reads]
    for s in strands:
        strand_hist[s] += 1


print("Total: %d reads, %d loci" % (total_reads, total_refstrs))
print("Passing filters: %d reads, %d loci" % (n_suppreads, n_refstrs))
print("Removed by X only filter: %d reads, %d loci" % (reads_x_only_filter, loci_x_only_filter))
print("Removed by Y only filter: %d reads, %d loci" % (reads_y_only_filter, loci_y_only_filter))
print("Removed by mapq60 filter: %d reads, %d loci" % (reads_mapq60_filter, loci_mapq60_filter))
print("Removed by min units filter: %d reads, %d loci" % (reads_min_units_filter, loci_min_units_filter))
print_sorted("Distribution of supp reads by chromosome", chrom_hist)
print_sorted("Distribution of supp reads by reference STR len", reflen_hist)
print_sorted("Distribution of supp reads by difference wrt reference STR len", reflen_diff_hist)
print_sorted("Distribution of reference STR loci by supp reads", nsupp_hist)
print_sorted("Distribution of reference STR loci by num of alleles", nalleles_hist)
print_sorted("Strand frequency over supp reads", strand_hist)
print_sorted("Mapping quality distribution over supp reads", mapq_hist)
