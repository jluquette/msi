#!/usr/bin/env python

import sys
import gzip
import time

if len(sys.argv) != 4:
    print('usage: %s sampleID mateID fastq.gz')
    exit(1)


class IlluminaFASTQIterator():
    """Reads FASTQs produced by Illumina CASAVA v1.8+ and splits them into
    chunks and read groups."""
    def __init__(self, filename, chunksize=4000000):
        if filename.lower().endswith('.gz'):
            self.f = gzip.open(filename, 'rb')
        else:
            self.f = open(filename, 'r')
        self.reads = 0
        self.chunk = 1
        self.chunksize = chunksize

    def __iter__(self):
        return self

    # FASTQ records are a 4 line series:
    #    @read name
    #    sequence
    #    +
    #    base quality scores
    def next(self):
        name = self.f.readline()
        seq = self.f.readline()
        dummy = self.f.readline()
        quals = self.f.readline()
        self.reads += 1

        if self.reads >= self.chunksize:
            self.reads = 0
            self.chunk += 1

        return (self.chunk, name, seq, dummy, quals)
        

# Illumina CASAVA 1.8+ read name format:
#    @HSQ700642:208:D1D6WACXX:2:1101:1199:2160 2:N:0:TGACCA
#    machineID:runID:flowcellID:laneID:tileID:xcoord:ycoord
#    mateID:fails filter:multiplex sequence
# We assume non-interleaved FASTQs (e.g., only one mateID will occur in
# the FASTQ and there will be at least two FASTQs per sample).
# We want files named:
#    sampleID_readgroup_mate_chunk.fastq
# Flowcell names are unique.  Lane+flowcell is a good read group definition.
def compute_readgroup(readname):
    main, aux = [ x.strip() for x in readname.split(' ') ]
    mainfields = [ x.strip() for x in main.split(':') ]
    return mainfields[2] + '-L00' + mainfields[3]


readgroups = {}
f = IlluminaFASTQIterator(sys.argv[3])
sample_id = sys.argv[1]
mate_id = sys.argv[2]
nreads = 0
clk = time.clock()
for chunk, name, seq, dummy, quals in f:
    readgroup = compute_readgroup(name)
    filename = "%s_%s_R%s_C00%d.fastq.gz" % (sample_id, readgroup, mate_id, chunk)
    try:
        outf = readgroups[filename]
    except KeyError:
        print("Opening new output file %s" % filename)
        outf = gzip.open(filename, 'wb')
        readgroups[filename] = outf

    outf.write(name)
    outf.write(seq)
    outf.write(dummy)
    outf.write(quals)

    nreads += 1
    if nreads % 1000000 == 0:
        print('wrote %d reads to %d files, %d seconds' % \
              (nreads, len(readgroups), time.clock() - clk))
        clk = time.clock()
