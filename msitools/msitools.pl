#!/usr/bin/env perl

# The previous version of this program split its original input into one
# file per chromosome and worked on a single chromosome at a time.  It also
# required the input to be sorted for efficient lookup of repeat records.
# This version is almost completely rewritten by Joe to neither split the
# input nor to require the input to be sorted.

# Using the ~8M repeat locus database produced by Tae-min, this program
# requires ~11GB of RAM to store the full database and human reference
# sequences in memory.  Typical RAM usage for analyzing the data in this
# experiment is ~50G.  One sample required as much as 52G.  If RAM is an
# issue, add option to specify a subset of the full chrlist and then
# merge the results in a separate script.

# TODO: TEST! compare against previous output.  will not be identical
# due to fixing the negative flanking index bug as well as removing all
# lines for repeats that have no supporting reads.

use strict;
use warnings;
use Getopt::Long;
use IO::Compress::Gzip qw($GzipError);
use lib 'perllib';
use Set::IntervalTree;
use Cwd;
use List::Util qw(max);


my $flank_bp = 2;
my $short_test = 0;
my $debug = 0;
my $resource_path = getcwd;
my @chrarray = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
undef my $inputfile;
undef my $outprefix;
GetOptions("flank_bp=i" => \$flank_bp,
           "input=s" => \$inputfile,
           "outprefix=s" => \$outprefix,
           "resource_path=s" => \$resource_path,
           "chr=s" => sub { @chrarray = map { s/chr//; $_; } split(",", $_[1]) },
           "debug" => \$debug)  # EXTREMELY verbose.  Do not use on large input
    or die("error parsing arguments");

die "--input is required" if not defined $inputfile;
die "--outprefix is required" if not defined $outprefix;
die "--resource_path is required" if not defined $resource_path;

print "flank_bp=$flank_bp\n";
print "input=$inputfile\n";
print "output prefix=$outprefix\n";

# Build a hash of accepted chromosomes
print "chroms=" . join(", ", @chrarray) . "\n";
my %chrhash;
foreach my $chr (@chrarray) {
    $chrhash{$chr} = 1;
}
print "Using " . scalar(keys %chrhash) . " chromosomes...\n";



# Reads the full set of human reference sequences into memory
print "Reading chromosome sequences into RAM..\n";
my %chrSeq;
foreach my $chr (keys %chrhash) {
    print "reading chromosome $chr.. ";
    open (F, "<", "$resource_path/chr$chr.fa") or die;
    my $null = <F>;  # strip off the first >chrXXX line
    my $nbytes = 0;
    while (<F>) {
        my $line = $_;
        chomp $line;
        $chrSeq{$chr} .= $line;
        $nbytes += length($line);
    }
    print "$nbytes bytes.\n";
    close F; 
}
print "done.\n";



# Build a record of all known repeats on this chromosome
# This used to be done one chrom at a time.  Changed by Joe.
print "Building repeat index..\n";
my %repeatdb;
my %repeatdb_stats;
open (F, "<", "$resource_path/WGRef_7892585MS_withGENEcategory_FINAL_withMSTYPE_hg19.txt") or die;
while (<F>) {
    chomp;
    my @element = split("\t");
    my $chr = $element[0];
    $chr =~ s/chr//g;
    # Only load repeats in the specified chromosome list
    if (exists $chrhash{$chr}) {
        $element[5] =~ s/\s//g;
        # There was a significant performance improvement for using an array
        # instead of a hash for these records, albeit at a cost to code
        # understandability.
        my $rec = [
            $chr,         # chr
            $element[1],  # start
            $element[2],  # end
            $element[3],  # region (exonic, intronic, intergenic, 3utr, 5utr)
            $element[4],  # seq
            $element[5],  # unit (mono, di, tri, tetra)
            uc(substr($chrSeq{$chr}, $element[1] - $flank_bp - 1, $flank_bp)),
            uc(substr($chrSeq{$chr}, $element[2] , $flank_bp)),  # flank2
            [],           # lens
            [],           # strands
            []            # mapqs
        ];
        # Intervals are half open: [a, b).
        $repeatdb{$chr} = Set::IntervalTree->new if not defined $repeatdb{$chr};
        $repeatdb{$chr}->insert($rec, $element[1], $element[2]+1);

        # The library doesn't make these values easily queryable.
        # track the number of records and the max record high value.
        $repeatdb_stats{$chr} = [ 0, 0 ] if not defined $repeatdb_stats{$chr};
        ++$repeatdb_stats{$chr}[0];
        $repeatdb_stats{$chr}[1] = max($repeatdb_stats{$chr}[1], $element[2]+1);
    }
}
close F;

while (my ($k, $v) = each(%repeatdb)) {
    print "chr$k: $repeatdb_stats{$k}[0] repeat records in [0, $repeatdb_stats{$k}[1])\n";
}
print "done.\n";



my $gzsuppreads = new IO::Compress::Gzip "$outprefix.supporting_reads.txt.gz"
    or die("IO::Compress::Gzip failed: $GzipError");
print $gzsuppreads "chr\tstart\tend\tregion\treadinfo\n";

my $nread = 0;
my $nsupp = 0;
my $num_searches = 0;
my $num_notindb = 0;
print "Reading input data..\n";
open (F, "<", $inputfile) or die;
while (<F>) {
    ++$nread;
    if ($nread % 10000 == 0) {
        print "$nread lines processed, $nsupp supporting reads, ";
        print $num_searches / 10000 . " mean searches per read.\n";
        $num_searches = 0;
    }

    chomp;
    my @element = split(" "); 
    my $chrom = $element[2];
    $chrom =~ s/chr//g;
    next if not exists $chrhash{$chrom};  # Skip if not in chr list

    my $startRepeat = $element[3] + $element[10];
    my $endRepeat = $element[3] + $element[11]; 
    my $strand = (int $element[1] & 0x10) == 0 ? '+' : '-';
    my $mapq = int $element[4];
    my $this_flank1 = substr($element[13], $element[10]-$flank_bp-1, $flank_bp);
    my $this_flank2 = substr($element[13], $element[11], $flank_bp);
    
    my $overlapping_repeats = $repeatdb{$chrom}->fetch($startRepeat, $endRepeat);
    if (scalar(@$overlapping_repeats) == 0) {
        if ($debug) {
            print "DEBUG: not in db: $chrom:[$startRepeat, $endRepeat)\n";
        }
        ++$num_notindb;
    }
    for my $record (@$overlapping_repeats) {
        my ($chr, $start, $end, $region, $seq, $unit, $flank1, $flank2, $lens, $strands, $mapqs) = @$record;

        ++$num_searches;
        if ($endRepeat >= $start and $startRepeat <= $end) {
            # joe: fix the case where $element[9]-$flank_bp-1 is negative,
            # which causes substr to select the end of the read
            if ($element[9] eq $unit and $flank1 eq $this_flank1 and
                $flank2 eq $this_flank2 and $element[10] - $flank_bp - 1 >= 0) {
                my $replen = scalar($element[11] - $element[10] + 1);
                push @$lens, $replen;
                push @$strands, $strand;
                push @$mapqs, $mapq;
                print $gzsuppreads "$chrom\t$start\t$end\t$region\t$replen\t@element\n";
                ++$nsupp;
                if ($debug) {
                    my $seq_context = uc(substr($chrSeq{$chrom},
                                                $start - $flank_bp - 1,
                                                $end - $start + 2*$flank_bp));
                    print "read: $chrom:$element[3], STR region: $startRepeat-$endRepeat\n";
                    print "record: $start\t$end\t$replen\t@element\n";
                    print "$flank1\t$flank2\n";
                    print "$this_flank1\t$this_flank2\n";
                    print "repeatSeq=" . " " x $flank_bp . "$seq\n";
                    print "chrSeq=   $seq_context\n";
                    print "lens here: @$lens\n";
                    print "strands here: @$strands\n";
                    print "mapqs here: @$mapqs\n";
                    print "-" x 80 . "\n";
                }
            }
            last; # Don't allow a read to match multiple repeat records
        }
    }
}
close F;
$gzsuppreads->close();
print "done.\n";
print "$nread total reads, $nsupp placed, $num_notindb not found in repeat db.\n";



print "Writing results to $outprefix.str_summary.txt.. ";
open (OUTPUT, ">$outprefix.str_summary.txt");
print OUTPUT "chr\tstart\tend\tunit\tregion\trepArray\tstrandArray\tmapQArray\n";
# Don't loop over keys %repeatdb, this way preserves chrom order
for my $chr (@chrarray) {  
    # Hack.  The IntervalTree library does not export a method for iterating
    # over its nodes, so we do a search that returns all possible nodes.
    # Recall that repeatdb_stats[1] is the max upper bound.
    my $all_recs = $repeatdb{$chr}->fetch(0, $repeatdb_stats{$chr}[1]);

    # Sort output records by start and end position (all recs have the same chr)
    for my $record (sort { $a->[1] <=> $b->[1] or $a->[2] <=> $b->[2] } @$all_recs) {
        my ($chr, $start, $end, $region, $seq, $unit, $flank1, $flank2, $lens, $strands, $mapqs) = @$record;
        print OUTPUT "$chr\t$start\t$end\t$unit\t$region";
        print OUTPUT "\t" . join(",", @$lens);
        print OUTPUT "\t" . join(",", @$strands);
        print OUTPUT "\t" . join(",", @$mapqs) . "\n";
    }
}
close OUTPUT;
print "done.\n";
