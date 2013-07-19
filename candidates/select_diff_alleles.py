#!/usr/bin/env python

"""
select interesting allele configurations:

We have 7 samples:
  * 2 bulk DNA from heart and cortex, PCR amplified
  * 1 100 neuron batch, MDA amplified
  * 4 single neurons, MDA amplified

"Interesting" allele configurations satisfy the following criteria:
  1. The 4 single neurons are either segregated by the allele or the
     entire set of single neurons is segregated from the heart tissue.
  2. The segregating allele appears EITHER:
     a. in heart bulk tissue alone, with no evidence in cortex bulk or
        the 100 neuron batch,
     b. in BOTH the cortex bulk and 100 neuron batch but not in the
        heart bulk tissue.
     By requiring the allele to be either present or absent in BOTH
     PCR and MDA amplified cortex samples, we hope to control for
     artifacts specific to MDA.
"""

import sys
from operator import itemgetter


if len(sys.argv) != 2:
    print('usage: %s candidate_file' % sys.argv[0])
    exit(1)

# Format is: allele1:depth,fracforward,meanmapq ... alleleN:... ...
# Return the set of alleles for which there is *any* evidence
def parse_alleles(summary):
    return set([ x.split(':')[0] for x in summary.split(' ') ])

def parse_genotype(gt):
    return set(gt.split('/'))

# Format is:
#   chrom, start, end, ref_len, unit, region
# followed by a variable number of 4-column groups:
#   #raw_alleles, call, genotype, allele_summary
f = open(sys.argv[1], 'r')
header_line = f.readline()
sys.stdout.write(header_line)
headers = [ h.strip() for h in header_line.split('\t') ]
heart_call_idx = headers.index('genotype.heart_bulk')
cortex_call_idx = headers.index('genotype.cortex_bulk')
neurons100_call_idx = headers.index('genotype.neurons_100batch')
neuron_call_idxs = \
    [ headers.index('genotype.' + x)
      for x in [ 'neuron_2', 'neuron_3', 'neuron_51', 'neuron_6' ] ]
heart_evidence_idx = headers.index('allele_summaries.heart_bulk')
cortex_evidence_idx = headers.index('allele_summaries.cortex_bulk')
neurons100_evidence_idx = headers.index('allele_summaries.neurons_100batch')
neuron_evidence_idxs = \
    [ headers.index('allele_summaries.' + x)
      for x in [ 'neuron_2', 'neuron_3', 'neuron_51', 'neuron_6' ] ]

for line in f:
    fields = [ x.strip() for x in line.split('\t') ]
    heart_alleles = parse_alleles(fields[heart_evidence_idx])
    cortex_alleles = parse_alleles(fields[cortex_evidence_idx])
    neurons100_alleles = parse_alleles(fields[neurons100_evidence_idx])
    neuron_alleles = parse_alleles(" ".join(itemgetter(*neuron_evidence_idxs)(fields)))

    any_brain = cortex_alleles | neurons100_alleles | neuron_alleles

    # Alleles which have supporting evidence in both PCR and MDA samples.
    # Can be from any of the MDA samples to alleviate dropout scenarios.
    brain_pcr_and_mda = cortex_alleles & (neurons100_alleles | neuron_alleles)

    # Step 1: find a segregating allele:
    # Do not consider any allele that has ANY evidence in either PCR or
    # MDA amplified brain tissue.
    seg_alleles = heart_alleles.symmetric_difference(any_brain)

    # By intersecting the PCR and MDA brain calls, we hope to remove alleles
    # that arose only under one amplification.  This could control for PCR
    # stutter alleles as well as allele dropout in MDA.
    seg_alleles &= brain_pcr_and_mda

    #print('heart: ' + str(heart_alleles))
    #print('cortex: ' + str(cortex_alleles))
    #print('100neuron: ' + str(neurons100_alleles))
    #print('single neurons: ' + str(neuron_alleles))
    #print('segregating: ' + str(alleles))
    
    if len(seg_alleles) == 0:
        continue

    # Now that we have a set of alleles that have some evidence, make sure
    # they were actually segregate in the genotype calls, which contain a
    # subset of the supported alleles.
    heart_calls = parse_genotype(fields[heart_call_idx])
    cortex_calls = parse_genotype(fields[cortex_call_idx])
    neurons100_calls = parse_genotype(fields[neurons100_call_idx])
    neuron_calls = reduce(lambda x, y: x | y,
                          [ parse_genotype(gt)
                            for gt in itemgetter(*neuron_call_idxs)(fields) ])
    print(heart_calls)
    print(cortex_calls)
    print(neurons100_calls)
    print(neuron_calls)
    
    sys.stdout.write(line)
    exit(1)
