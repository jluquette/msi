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
  * STR unit distribution: mono, di, tri, tetra
  * STR unit multiple distribution: (end-start) / unit
  * STR genomic location distribution: intergenic, intronic, exonic

Options for filtering reads:
  * haploid options: --only-x, --only-y, only use X or Y alignments, depending
    on option (both cannot be specified together).  For male samples, these
    should be haploid and thus deviations from expectation should reflect
    experimental errors.
  * maximum mapQ option: --mapq60, only use mapq60 (very high confidence)
    alignments
  * minimum repeat unit option: --min-unit (recommended=3, default=1?).  Only
    consider loci where the repeat unit (end-start+1) / unit size is greater
    than the specified value.  This could be useful to remove questionable
    loci, like 2~3 repeat units of tri or tetranucleotide repeats.
"""

import sys
from argparse import ArgumentParser
from strlocusiterator import STRLocusIterator

parser = ArgumentParser()
parser.add_argument('filename', metavar='str_summary', type=str,
                    help='STR summary file from msitools')
parser.add_argument('--min-mapq', dest='min_mapq', metavar='N', default=0,
                    type=int,
                    help='Discard reads with mapping quality < N')
parser.add_argument('--min-units', dest='min_units', metavar='N', default=1,
                    type=int,
                    help='Discard reference loci with < N repeat units')
parser.add_argument('--min-supp', dest='min_supp_reads', metavar='N', default=0,
                    type=int,
                    help='Discard reference loci with < N supporting reads ' \
                         'after all read filters have been applied')
parser.add_argument('--max-ref-diff', dest='max_ref_diff', metavar='N',
                    default=0, type=int,
                    help='Discard reads that differ too greatly from the ' \
                         'reference STR length.  abs(observed len - ref len) ' \
                         '> N, for N in base pairs')
parser.add_argument('--x-only', dest='x_only', action='store_true',
                    default=False,
                    help='Only consider reads from the X chromosome')
parser.add_argument('--y-only', dest='y_only', action='store_true',
                    default=False,
                    help='Only consider reads from the Y chromosome')

args = parser.parse_args()

if args.x_only and args.y_only:
    print('ERROR: only one of --x-only and --y-only may be specified')
    exit(1)

locus_f = STRLocusIterator(**vars(args))
for (chrom, start, end, unit, region, reads) in locus_f:
    # Don't do anything, just accumulate metrics
    continue

for (description, value) in locus_f.filter_metrics():
    print("%s\t%d" % (description, value))

for (description, hist) in locus_f.hist_metrics():
    print(description)
    for k in sorted(hist.keys()):
        print("%s\t%d" % (k, hist[k]))
