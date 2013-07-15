MSITools is a set of perl scripts developed by Tae-min for mapping Sputnik's
results to a library of known repeat regions in the genome and for summarizing
the mapping.

Changes to Tae-min's code:
    * fix left-flank loop-to-beginning of read bug
    * change the flank_size requirement from 2bp->10bp
    * msitools now creates a file of all reads that were used to support one
      of its calls in a file named 'SUPPORT_reads_<filename>.txt'
