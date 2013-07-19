#!/usr/bin/env python

"""
Print out loci and read summaries for sites passing several filter criteria.
Differs from metrics_msitools primarily in that it outputs the loci, and in
the future will compute some additional information like a genotype call with
some confidence metric.
"""

import sys
from argparse import ArgumentParser
from strlocusiterator import STRLocusIterator
from collections import defaultdict
from numpy import array

parser = ArgumentParser()
STRLocusIterator.add_parser_args(parser)
args = parser.parse_args()

def summarize_reads(reads):
    """`reads` is a list of (obs len, strand, mapq) tuples each describing
    a single read.  Return 3 summary strings, one for observed lengths, one
    summarizing strands and one summarizing mapping qualities."""
    # Map: observed len -> (num reads, num forward strand, sum mapq)
    metrics = defaultdict(lambda: array([0, 0, 0]))
    for obs, strand, mapq in reads:
        metrics[obs] += array([ 1, 1 if strand == '+' else 0, mapq ])

    return " ".join("%d:%d,%.2f,%.2f" % \
                    (k, n, float(nforward) / n, float(summapq) / n)
                        for k, (n, nforward, summapq) in metrics.iteritems())


locus_f = STRLocusIterator(**vars(args))
for (chrom, start, end, unit, region, reads) in locus_f:
    # Don't do anything, just accumulate metrics
    allele_summary = summarize_reads(reads)
    print "\t".join([ chrom, str(start), str(end), unit, region, allele_summary])
