Indexed the same FASTAs used to construct Tae-min's repeat database.
Concatenated the 24 files into a single hg19.fa.  No decoy sequence, no
unplaced contigs, no hapmap alternate haplotypes, etc..

../msitools/hg19_chromosome_seq/chr{1..22,X,Y}.fa

../dependencies/bwa-0.7.5a/bwa index

Another modification to Tae-min's WGRef database: the following line has
no sequence associated with it, so it was removed.

1691765 chr3 144348257 144348268 intergenic                mono

Here's a curious entry in Tae-min's database:

Y       2762709 2762721 13      mono    intergenic      TTCTTTCTTT      AAGACAGGGT      TTTTTTTTTTTTT

Notice that the repeat is (T)n and that its left flank has a string of 3 Ts
on its right border.  Why aren't those 3 Ts included in the repeat sequence?
I checked on UCSC genome browser, and this is indeed the correct sequence.
