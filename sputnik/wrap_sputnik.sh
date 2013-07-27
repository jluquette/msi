#!/bin/bash

# This script wraps Sputnik to give the illusion that it can both accept
# a BAM file as input and produce a BAM file as output.  The final BAM
# file will include special tags in the optional SAM tag positions to
# encode Sputnik's repeat sensing.  Reads which do not contain STRs will
# be not be included in the final BAM.

if [ $# != 2 ]; then
    echo "usage: $0 input_bam output_bam"
    exit 1
fi

input=$1
output=$2

# The header goes out unmodified
samtools view -H $input > /dev/null

# Sputnik expects a FASTA style input file.  Input is two lines per read,
# one description line starting with > and then a line containing the read
# sequence.  We piggyback the entire SAM line onto the description line and
# then extract only the sequence for the second line.  Example:
# SAM input line:
#    HSQ700642:191:D18JJACXX:5:1103:4330:4941       83      chr1    28517   24     74M3I23M =       28505   -109    TTTCCTGCCTATCCATTTTGTTAACTCTTCAATGCATTCCACAAATGCCTAAGTATTCTTTAATAATGGTGGTTTGGTTTTTTTTTTTTTGCATCTATGA    ####?@CCCC>::3:BDCC>C>44+(3>>3C@>@>@A@<9;>CCC>@DC@>@;;;A@9CEBC@;BBB?;6;@@EEAGGJJJJJJJJJHGHHHFFFFFCCC   RG:Z:D18JJACXX-L005      NM:i:4  AS:i:83 XS:i:80
# Translates via awk to:
#    >HSQ700642:191:D18JJACXX:5:1103:4330:4941       83      chr1    28517   24     74M3I23M =       28505   -109    TTTCCTGCCTATCCATTTTGTTAACTCTTCAATGCATTCCACAAATGCCTAAGTATTCTTTAATAATGGTGGTTTGGTTTTTTTTTTTTTGCATCTATGA    ####?@CCCC>::3:BDCC>C>44+(3>>3C@>@>@A@<9;>CCC>@DC@>@;;;A@9CEBC@;BBB?;6;@@EEAGGJJJJJJJJJHGHHHFFFFFCCC   RG:Z:D18JJACXX-L005      NM:i:4  AS:i:83 XS:i:80
#    TTTCCTGCCTATCCATTTTGTTAACTCTTCAATGCATTCCACAAATGCCTAAGTATTCTTTAATAATGGTGGTTTGGTTTTTTTTTTTTTGCATCTATGA
# Sputnik's output is:
#    description line {mono|di|tri|tetra} repeatStart repeatEnd score sequence
# e.g.:
#    >HSQ700642:191:D18JJACXX:5:1103:4330:4941       83      chr1    28517   24     74M3I23M =       28505   -109    TTTCCTGCCTATCCATTTTGTTAACTCTTCAATGCATTCCACAAATGCCTAAGTATTCTTTAATAATGGTGGTTTGGTTTTTTTTTTTTTGCATCTATGA    ####?@CCCC>::3:BDCC>C>44+(3>>3C@>@>@A@<9;>CCC>@DC@>@;;;A@9CEBC@;BBB?;6;@@EEAGGJJJJJJJJJHGHHHFFFFFCCC   RG:Z:D18JJACXX-L005      NM:i:4  AS:i:83 XS:i:80 mono 78 90 12 TTTCCTGCCTATCCATTTTGTTAACTCTTCAATGCATTCCACAAATGCCTAAGTATTCTTTAATAATGGTGGTTTGGTTTTTTTTTTTTTGCATCTATGA


# Since perl isn't my strong suit either: -a enables autosplit mode, which
# splits the input line by the -F delimiter.  -n wraps the statement in a
# while (<>) { ... } loop.
# IMPORTANT!!! Sputnik's 5 appended fields are SPACE separated, not tab
# separated.  This means the final field (by tab separation) will be the
# last optional SAM tag (is there always at least one tag?) plus all of
# Sputnik's fields.
# The inline perl program:
#   * removes the ">" string from the read name, which is added by awk in
#     the FASTA conversion step,
#   * prints the first through (last-1)th SAM fields
#   * recognizes that the final SAM field contains BOTH a SAM tag AND all
#     of Sputnik's added annotations; all of these elements are separated
#     by spaces.  It then splits the final field by space, outputs the final
#     SAM tag and then outputs 4 special tags for storing the Sputnik data:
#     ru=repeat unit, rs=repeat start, re=repeat end, ss=Sputnik score.
#     (The SAM v1.4 spec says that all lower case tag names are reserved
#     for end users.)
(samtools view -H $input; samtools view $input \
    | awk '{ print ">" $0; print $10; }' \
    | ./sputnik_V6_len5 /dev/stdin \
    | perl -lan -F/\\t/ -e 'my @sput=split(" ", $F[$#F]); print substr($F[0], 1) . "\t" . join("\t", @F[1 .. ($#F - 1)]) . "\t$sput[0]\tru:Z:$sput[1]\trs:i:$sput[2]\tre:i:$sput[3]\tss:i:$sput[4]"') \
    | samtools view -bSo $output -