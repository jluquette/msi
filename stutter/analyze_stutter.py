#!/usr/bin/env python

import sys
from scipy.stats import binom_test

if len(sys.argv) != 2:
    print("usage: %s candidates.txt" % sys.argv[0])
    exit(1)

binom_threshold = 0.05
total = 0
reads_total = 0
n_homozygous = 0
reads_homozygous = 0
n_probably_stutter = 0
reads_probably_stutter = 0
n_probably_heterozygous = 0
reads_probably_heterozygous = 0
with open(sys.argv[1], 'r') as f:
    for line in f:
        fields = [ x.strip() for x in line.split('\t') ]
        # allele len:depth,frac forward,mean mapq
        allele_info = [ x.split(':') for x in fields[10].split(' ') ]
        # (len, depth, frac forward, mean mapq)
        alleles = [ (x[0],) + tuple(x[1].split(',')) for x in allele_info ]

        total += 1
        for x in alleles:
            reads_total += int(x[1])

        # Use the "best" two alleles (best=highest read count) and assume
        # equal probability of occurence.  Then model with a binomial and
        # get a p-value. (we use this "best two" approximation since scipy
        # does not have a built in multinomial model)
        if len(alleles) > 1:
            # sort by index 1 in the tuple (read depth) and take the top 2
            best_alleles = sorted(alleles, key=lambda x: int(x[1]), reverse=True)[0:2]
            # grab the depth for each of the best alleles
            #print(best_alleles)
            obs_a, obs_b = [ int(x[1]) for x in best_alleles ]
            #print("a=%d, b=%d" % (obs_a, obs_b))
            pval = binom_test(obs_a, n=obs_a + obs_b, p=0.5)
            if pval < binom_threshold:
                n_probably_stutter += 1
                # should really be all obs reads != max(obs_a, obs_b)
                reads_probably_stutter += min(obs_a, obs_b)
                print("stutter\t%0.5f\t" % pval + line.strip())
            else:
                reads_probably_heterozygous += obs_a + obs_b
                n_probably_heterozygous += 1
                print("het\t%0.5f\t" % pval + line.strip())
        else:
            n_homozygous += 1
            for x in alleles:
                reads_homozygous += int(x[1])


print("binom threshold: %.3f" % binom_threshold)
print("hom: %d" % n_homozygous)
print("het: %d" % n_probably_heterozygous)
print("stutter: %d" % n_probably_stutter)
print("total: %d" % total)
print("reads hom: %d" % reads_homozygous)
print("reads het: %d" % reads_probably_heterozygous)
print("reads stutter: %d" % reads_probably_stutter)
print("reads total: %d" % reads_total)
