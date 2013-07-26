Indexed the same FASTAs used to construct Tae-min's repeat database.
Concatenated the 24 files into a single hg19.fa.  No decoy sequence, no
unplaced contigs, no hapmap alternate haplotypes, etc..

../msitools/hg19_chromosome_seq/chr{1..22,X,Y}.fa

../dependencies/bwa-0.7.5a/bwa index
