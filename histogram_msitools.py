#!/usr/bin/env python

"""
histogram_sputnik.py - Make a frequency table of #supporting reads per locus.
"""

import sys
from collections import defaultdict

# line format is tab-delimited: index, chrom, start, end, read count list
# read count list is a string of numbers separated by commas.  The numbers
# are the repeat unit lengths; each number represents a single read.  So a
# string like 7,7,7,8 means 3 reads support an allele with 7 repeat units
# and a single read supports an allele with 8 repeat units.
# We drop the 'index' field.
freq_counts = defaultdict(int)
f = open(sys.argv[1], 'r')
for line in f:
    fields = map(str.strip, line.split('\t'))

    # header row
    if fields[0] == 'index':
        continue

    if not fields[4]:
        # if the read count list is empty, there are no reads supporting
        # any repeat length allele.  skip these.
        continue

    # filter(None, seq) drops false values in seq.  Empty strings in
    # python are "falsey": http://docs.python.org/2/library/stdtypes.html#truth-value-testing
    supports = filter(None, map(str.strip, fields[4].split(',')))
    freq_counts[len(supports)] += 1

for k,v in freq_counts.items():
    print(str(k) + '\t' + str(v))
