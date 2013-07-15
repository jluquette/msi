#!/usr/bin/env python

"""Insert the lobSTR header between #-started lines and data lines."""

import sys

if len(sys.argv) != 2:
    print('usage: %s lobstr_genotype_tab' % sys.argv[0])
    exit(1)

header_f = open('header.txt', 'r')
header = header_f.readline().split()
header_f.close()

lobstr_f = open(sys.argv[1].strip(), 'r')
in_comments = True
for line in lobstr_f:
    if in_comments and line.strip()[0] != '#':
        in_comments = False
        print('\t'.join(header))

    if line.strip()[0] == '#':
        print(line),
    else:
        print('\t'.join(line.split()))
