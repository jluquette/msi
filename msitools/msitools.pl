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
        # There was a significant performance improvement for using an array
        # instead of a hash for these records, albeit at a cost to code
        # understandability.
        my $rec = [
            $chr,         # chr
            $element[1],  # start
            $element[2],  # end
            $element[4],  # seq
            $element[5],  # type
            uc(substr($chrSeq{$chr}, $element[1] - $flank_bp - 1, $flank_bp)),
            uc(substr($chrSeq{$chr}, $element[2] , $flank_bp)),  # flank2
            [],           # lens
            [],           # strands
            []            # mapqs
        ];
        push @{ $repeatdb{$chr} }, $rec;
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
my $num_searches = 0;
my $idx;
my $chrom;
print "Reading input data..\n";
open (F, "<", $inputfile) or die;
while (<F>) {
    ++$nread;
    if ($nread % 10000 == 0) {
        print "$nread lines processed, $nsupp supporting reads, ";
        print $num_searches / $nread . " mean searches per read.\n";
        $num_searches = 0;
    }

    chomp;
    my @element = split(" "); 
    $element[2] =~ s/chr//g;

    # Optimize the search: if $chrom hasn't changed, then don't reset $idx.
    # This will significantly increase performance when input is sorted.
    $idx = 0 if ($chrom ne $element[2]);
    $chrom = $element[2];
    next if not exists $chrhash{$chrom};  # Skip if not recognized

    my $rec_array = $repeatdb{$chrom};
    my $direction = ""; # searching direction
    my $startRepeat = $element[3] + $element[10];
    my $endRepeat = $element[3] + $element[11]; 
    my $strand = (int $element[1] & 0x10) == 0 ? '+' : '-';
    my $mapq = int $element[4];
    my $this_flank1 = substr($element[13], $element[10]-$flank_bp-1, $flank_bp);
    my $this_flank2 = substr($element[13], $element[11], $flank_bp);
    while (1) { # search all repeats for an overlapping record
        #my $record = $rec_array->[$idx];  # ref to a hash
        my ($chr, $start, $end, $seq, $type, $flank1, $flank2, $lens, $strands, $mapqs) = @{ $rec_array->[$idx] };
        if ($endRepeat >= $start and $startRepeat <= $end) {
            # joe: fix the case where $element[9]-$flank_bp-1 is negative,
            # which causes substr to select the end of the read
            if ($element[9] eq $type and $flank1 eq $this_flank1 and
                $flank2 eq $this_flank2 and $element[10] - $flank_bp - 1 >= 0) {
                my $replen = scalar($element[11] - $element[10] + 1);
                push @$lens, $replen;
                push @$strands, $strand;
                push @$mapqs, $mapq;
                print $gzsuppreads "$chrom\t$start\t$end\t$replen\t@element\n";
                ++$nsupp;
                if ($debug) {
                    my $seq_context = uc(substr($chrSeq{$chrom},
                                                $start - $flank_bp - 1,
                                                $end - $start + 2*$flank_bp));
                    print "$start\t$end\t$replen\t@element\n";
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
        } elsif ($startRepeat > $end) { 
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
                ++$num_searches;
            } 
        } elsif ($endRepeat < $start) {
            # Backward searching (Read repeat >> Genome repeat)
            if ($direction eq "forward" or $idx <= 0) {
                last;
            } else {
                $direction = "backward";
                --$idx;
                ++$num_searches;
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
        my ($chr, $start, $end, $seq, $type, $flank1, $flank2, $lens, $strands, $mapqs) = @$record;
        print OUTPUT "$chr\t$start\t$end";
        print OUTPUT "\t" . join(",", @$lens);
        print OUTPUT "\t" . join(",", @$strands);
        print OUTPUT "\t" . join(",", @$mapqs) . "\n";
    }
}
close OUTPUT;
print "done.\n";
