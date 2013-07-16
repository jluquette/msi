#!/usr/bin/env perl

# TODO: TEST! compare against previous output.  will not be identical
# due to fixing the negative flanking index bug as well as removing all
# lines for repeats that have no supporting reads.

use strict;
use warnings;
use Getopt::Long;
use IO::Compress::Gzip qw($GzipError);

my $flank_bp = 2;
my $short_test = 0;
my $debug = 0;
undef my $inputfile;
undef my $outprefix;
undef my $resource_path;
GetOptions("flank_bp=i" => \$flank_bp,
           "input=s" => \$inputfile,
           "outprefix=s" => \$outprefix,
           "resource_path=s" => \$resource_path,
           "test" => \$short_test,
           "debug" => \$debug)  # EXTREMELY verbose.  Do not use on large input
    or die("error parsing arguments");

die "--input is required" if not defined $inputfile;
die "--outprefix is required" if not defined $outprefix;
die "--resource_path is required" if not defined $resource_path;

print "flank_bp=$flank_bp\n";
print "input=$inputfile\n";
print "output prefix=$outprefix\n";

# Build a hash of accepted chromosomes
my @chrarray = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
if ($short_test) {
    @chrarray = qw(1 Y);
}
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
    open (F, "<", "$resource_path/hg19_chromosome_seq/chr$chr.fa") or die;
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
open (F, "<", "$resource_path/WGRef_7892585MS_withGENEcategory_FINAL_withMSTYPE_hg19.txt") or die;
while (<F>) {
    chomp;
    my @element = split("\t");
    my $chr = $element[0];
    $chr =~ s/chr//g;  # Get rid of 'chr' if it's there
    # Only load repeats in the specified chromosome list
    if (exists $chrhash{$chr}) {
        $element[5] =~ s/\s//g;
        my %rec = (
            'chr'     => $chr,
            'start'   => $element[1],
            'end'     => $element[2],
            'seq'     => $element[4],
            'type'    => $element[5],
            'flank1'  => uc(substr($chrSeq{$chr},
                                   $element[1] - $flank_bp - 1,
                                   $flank_bp)),
            'flank2'  => uc(substr($chrSeq{$chr}, $element[2] , $flank_bp)),
            'lens'    => [],
            'strands' => [],
            'mapqs'   => []
        );
        push @{ $repeatdb{$chr} }, { %rec };
    }
}
close F;

while (my ($k, $v) = each(%repeatdb)) {
    print "chr$k: " . scalar(@$v) . " repeat records\n";
}
print "done.\n";



my $gzsuppreads = new IO::Compress::Gzip "$outprefix.supporting_reads.txt.gz"
    or die("IO::Compress::Gzip failed: $GzipError");
print $gzsuppreads "chr\tstart\tend\tSTRlen\treadinfo\n";

my $nread = 0;
my $nsupp = 0;
my $idx;
my $last_chrom = '';
print "Reading input data..\n";
open (F, "<", $inputfile) or die;
while (<F>) {
    ++$nread;
    if ($nread % 10000 == 0) {
        print "$nread lines processed, $nsupp supporting reads\n";
    }

    chomp;
    my @element = split(" "); 
    my $chrom = $element[2];
    $chrom =~ s/chr//g;                   # Get rid of 'chr' if it's there
    next if not exists $chrhash{$chrom};  # Skip if not recognized
    $idx = 0 if ($chrom ne $last_chrom);
    my $rec_array = $repeatdb{$chrom};
    my $direction = ""; # searching direction

    my $startRepeat = $element[3] + $element[10];
    my $endRepeat = $element[3] + $element[11]; 
    my $strand = (int $element[1] & 0x10) == 0 ? '+' : '-';
    my $mapq = int $element[4];
    my $this_flank1 = substr($element[13], $element[10]-$flank_bp-1, $flank_bp);
    my $this_flank2 = substr($element[13], $element[11], $flank_bp);
    while (1) { # search all repeats for an overlapping record
        my $record = $rec_array->[$idx];  # ref to a hash
        if ($endRepeat >= $record->{start} and $startRepeat <= $record->{end}) {
            # joe: fix the case where $element[9]-$flank_bp-1 is negative,
            # which causes substr to select the end of the read
            if ($element[9] eq $record->{type} and
                $record->{flank1} eq $this_flank1 and
                $record->{flank2} eq $this_flank2 and
                $element[10] - $flank_bp - 1 >= 0) {
                my $replen = scalar($element[11] - $element[10] + 1);
                push @{ $record->{lens} }, $replen;
                push @{ $record->{strands} }, $strand;
                push @{ $record->{mapqs} }, $mapq;
                print $gzsuppreads "$chrom\t$record->{start}\t$record->{end}\t$replen\t@element\n";
                ++$nsupp;
                if ($debug) {
                    my $seq_context = 
                        substr($chrSeq{$chrom},
                               $record->{start} - $flank_bp - 1,
                               $record->{end} - $record->{start} + 2*$flank_bp));
                    print "$record->{start}\t$record->{end}\t$replen\t@element\n";
                    print "$record->{flank1}\t$record->{flank2}\n";
                    print "$this_flank1\t$this_flank2\n";
                    print "repeatSeq=" . " " x $flank_bp . "$record->{seq}\n";
                    print "chrSeq=   " . uc($seq_context) . "\n"
                    print "lens here: @{ $record->{lens} }\n";
                    print "strands here: @{ $record->{strands} }\n";
                    print "mapqs here: @{ $record->{mapqs} }\n";
                    print "-" x 80 . "\n";
                }
            }
            last; # Don't allow a read to match multiple repeat records
        } elsif ($startRepeat > $record->{end}) { 
            # Forward searching (Genome repeat << Read repeat)
            # Fix by Joe: used to be $indexPos >= $repeatno, but this allows
            # ++indexPos to be run when indexPos = repeatno - 1, which means
            # indexPos will become repeatno, which is not a valid index to
            # repeatStart, repeatEnd, etc. (repeatno is now scalar(@$rec_array))
            if ($direction eq "backward" or $idx >= scalar(@$rec_array) - 1) {
                last;
            } else {
                $direction = "forward";
                ++$idx;
            } 
        } elsif ($endRepeat < $record->{start}) {
            # Backward searching (Read repeat >> Genome repeat)
            if ($direction eq "forward" or $idx <= 0) {
                last;
            } else {
                $direction = "backward";
                --$idx;
            }
        } else {
            print "ERROR: no repeat record found for read:\n";
            die("$chrom:$element[0]-$element[2]: $startRepeat-$endRepeat");
        }
    }
}
close F;
$gzsuppreads->close();
print "done.\n";

print "Writing results to $outprefix.str_summary.txt.. ";
open (OUTPUT, ">$outprefix.str_summary.txt");
print OUTPUT "index\tchr\tstart\tend\trepArray\tstrandArray\tmapQArray\n";
# Don't loop over keys %repeatdb, this way preserves chrom order
for my $chr (@chrarray) {  
    for my $record (@{ $repeatdb{$chr} }) {
        if (scalar(@{ $record->{lens} }) > 0) {
            print OUTPUT "$chr\t$record->{start}\t$record->{end}";
            print OUTPUT "\t" . join(",", @{ $record->{lens} });
            print OUTPUT "\t" . join(",", @{ $record->{strands} });
            print OUTPUT "\t" . join(",", @{ $record->{mapqs} }) . "\n";
        }
    }
}
close OUTPUT;
print "done.\n";
