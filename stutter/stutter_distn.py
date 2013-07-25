#!/usr/bin/env python

import sys
from numpy import array
from collections import defaultdict
from operator import itemgetter
from scipy.stats import binom_test

if len(sys.argv) != 2:
    print('usage: %s str_summary_file' % sys.argv[0])
    exit(1)

unit_to_int = { 'mono':1, 'di':2, 'tri':3, 'tetra':4 }

def parse_summary(s):
    lendiff, info = s.split(':')
    nreads, fracforward, meanmapq = info.split(',')
    return (int(lendiff), int(nreads), float(fracforward), float(meanmapq))


def get_stats(summaries, use_binom=False, binom_threshold=0.05):
    ainfo = [ parse_summary(s) for s in summaries.split(' ') ]
    # Sort by highest number of reads
    sorted_ainfo = sorted(ainfo, key=lambda x: x[1], reverse=True)

    likely_stutter = 0
    if use_binom:
        # Use a simple but arbitrary binomial model to distinguish between
        # heterozygous loci and polymerase stutter
        if len(sorted_ainfo) > 1:
            # Get the "best" two allele lens
            obs1 = sorted_ainfo[0][1]
            obs2 = sorted_ainfo[1][1]
            if binom_test(obs1, n=obs1 + obs2, p=0.5) <= binom_threshold:
                # Reject hypothesis that each allele is equally likely: stutter
                likely_stutter = obs2    # obs1 is never considered stutter
                #print("stutter: " + summaries)
    
        # reads supporting any alleles other than the top 2 are automatically
        # considered stutter
        likely_stutter += sum([ x[1] for x in sorted_ainfo[2:] ])
    else:
        # Assume EVERYTHING but the primary allele is due to polymerase
        # stutter.  NOTE: even on hemizygous sex chromosomes, this assumption
        # can be invalid for, e.g., bulk PCR samples that profile several
        # individual cells, each of which may have its own allele.
        likely_stutter = sum([ x[1] for x in sorted_ainfo[1:] ])

    tot_reads = sum([ x[1] for x in ainfo ])
    return array([ tot_reads, likely_stutter, 1, 1 if likely_stutter else 0 ])


# Values:
#    [ total reads, likely stutter reads, total loci, likely stutter loci ]
stats = defaultdict(lambda: defaultdict(lambda: array([0, 0, 0, 0])))

# File format:
#     chr, start, end, ref STR len (bp), unit (mono, di, tri, tetra), region,
#     #alleles before genotyping, call, genotype, pval, allele summaries
# The allele summary format: space separated list of:
#     allele len diff from reference:depth,frac forward,meanmapq
with open(sys.argv[1], 'r') as f:
    f.readline()  # skip header
    for line in f:
        fields = [ x.strip() for x in line.split('\t') ]
        reflen, unit, summaries = itemgetter(3,4,10)(fields)
        stats[unit][int(reflen)] += get_stats(summaries, use_binom=True)


print("unit\treflen\tstutter_reads\ttotal_reads\tpercent_stutter_reads\tstutter_loci\ttotal_loci\tpercent_stutter_loci")
for unit in sorted(stats.keys(), key=lambda x: unit_to_int[x]):
    for reflen in sorted(stats[unit].keys()):
        x = stats[unit][reflen]
        print("%s\t%d\t%d\t%d\t%.3f\t%d\t%d\t%.3f" % \
              (unit, reflen, x[1], x[0], float(x[1])/x[0],
               x[3], x[2], float(x[3])/x[2]))
