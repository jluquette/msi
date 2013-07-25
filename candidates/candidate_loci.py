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
from itertools import product
from numpy import array
from scipy.stats import binom
from numpy import prod

parser = ArgumentParser()
STRLocusIterator.add_parser_args(parser)
args = parser.parse_args()

# Only in python 2.7+...
def combinations_with_replacement(iterable, r):
    pool = tuple(iterable)
    n = len(pool)
    for indices in product(range(n), repeat=r):
        if sorted(indices) == list(indices):
            yield tuple(pool[i] for i in indices)


def summarize_alleles(reads, reflen):
    """`reads` is a list of (obs len, strand, mapq) tuples each describing
    a single read.
    `reflen` is the length of the STR in the reference genome."""
    # Map: observed len -> (num reads, num forward strand, sum mapq)
    metrics = defaultdict(lambda: array([0, 0, 0]))
    for obs, strand, mapq in reads:
        metrics[obs - reflen] += array([ 1, 1 if strand == '+' else 0, mapq ])

    return dict((k, (n, float(nforward) / n, float(summapq) / n))
                for k, (n, nforward, summapq) in metrics.iteritems())


def genotype(alleles, total_err=0.01):
    """`total_err` is totally arbitrary.  It is meant to model the probability
    of any kind of error producing an observed read.  It must account for
    alleleic imbalance, polymerase stutter, sequencing error, alignment error,
    etc.  Its current value is purely arbitrary."""

    # Try every possible selection of 2 alleles from the set
    gts = list(combinations_with_replacement(alleles.keys(), 2))

    def hap_prob(x, Nalleles, haptype, err):
        return 1 - err if x == haptype else err/(Nalleles - 1)

    def dip_prob(alleles, hapA, hapB, err):
        return prod([ (hap_prob(k, len(alleles), hapA, err)/2 + \
                       hap_prob(k, len(alleles), hapB, err)/2)**n
                      for k, (n, _1, _2) in alleles.iteritems() ])

    probs = [ dip_prob(alleles, a, b, total_err) for a, b in gts ]

    best = gts[probs.index(max(probs))]
    if best[0] == best[1]:
        if best[0] == 0:
            call = 'ref'
        else:
            call = 'hom'
    else:
        call = 'het'

    return (call, max(probs), str(best[0]) + "/" + str(best[1]))


locus_f = STRLocusIterator(**vars(args))
print("chr\tstart\tend\tref_len\tunit\tregion\traw_alleles\tcall\tgenotype\tpval\tallele_summaries")
for (chrom, start, end, unit, region, reads) in locus_f:
    # Don't do anything, just accumulate metrics
    reflen = end - start + 1
    allele_summaries = summarize_alleles(reads, reflen)
    allele_strings = " ".join("%d:%d,%.2f,%.2f" % ((k,) + v)
                              for k, v in allele_summaries.iteritems())
    call, pval, gt = genotype(allele_summaries)
    print("\t".join([ chrom, str(start), str(end), str(reflen), unit, region,
                      str(len(allele_summaries)), call, gt,
                      "%g" % pval, allele_strings ]))
