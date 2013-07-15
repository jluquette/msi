Original BAMs lived here.  Aligned with standard BWA aln (not mem), version
0.6.2-r126 by Gilad.  PCR duplicates were removed by Picard and pairs of
mapQ=0 reads were removed by sam_to_sputnik.py, which produces a FASTA file
from the BAM (after conversion to SAM via samtools).

The markdup.log and markdup_metrics.txt files show the effect of PCR duplicate
removal.  The .nlines.txt files show the effect of removing mapQ=0 pairs AND
converting to FASTA.

NOTE NOTE NOTE: since Sputnik takes FASTA as input, it does not use any base
quality score data to inform its calls.  It also means that the base quality
data is lost for downstream analysis unless we do extra work to add it back
in.
