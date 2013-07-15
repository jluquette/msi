#!/usr/bin/env perl

# TODO:
#  make flank_size a command line arg
#  use a better perl argument parsing library
#  fix the field indeces in the parsing at the end of the script.  (i.e.,
#  since I added the mapQ value to the line, everything needs to shift over
#  by one.)

use strict;
use warnings;
use Getopt::Long;

my $flank_bp = 2;
my $short_test = 0;
undef my $inputfile;
undef my $outprefix;
undef my $resource_path;
GetOptions("flank_bp:i" => \$flank_bp,
           "input=s" => \$inputfile,
           "outprefix=s" => \$outprefix,
           "resource_path=s" => \$resource_path,
           "test" => \$short_test)
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
        $chrSeq{"chr$chr"} .= $line;
        $nbytes += length($line);
    }
    print "$nbytes bytes.\n";
	close F; 
}
print "done.\n";



# Build a record of all known repeats on this chromosome
# This used to be done one chrom at a time.  Changed by Joe.
print "Building repeat index..\n";
my @repeatChrom = ();
my @repeatStart = ();
my @repeatEnd = ();
my @repeatSeq = ();
my @repeatFlank1 = ();
my @repeatFlank2 = ();
my @repeatType = ();
my $repeatno = 0;
open (F, "<", "$resource_path/WGRef_7892585MS_withGENEcategory_FINAL_withMSTYPE_hg19.txt") or die;
while (<F>) {
    my $line = $_;
	chomp $line;
    my @element = split("\t", $line); 
    my $chr = $element[0];
    $chr =~ s/chr//g;  # Get rid of 'chr' if it's there
    # Only load repeats in the specified chromosome list
    if (exists $chrhash{$chr}) {
	    $element[5] =~ s/\s//g;
        $repeatChrom[$repeatno] = $element[0];
        $repeatStart[$repeatno] = $element[1]; 
	    $repeatEnd[$repeatno] = $element[2];
        $repeatSeq[$repeatno] = $element[4];
	    $repeatType[$repeatno] = $element[5];
	    $repeatFlank1[$repeatno] =
            uc(substr($chrSeq{"chr".$chr}, $element[1]-$flank_bp-1, $flank_bp));
	    $repeatFlank2[$repeatno] =
            uc(substr($chrSeq{"chr".$chr}, $element[2] , $flank_bp));
	    ++$repeatno;
    }
}
close F;
print "done.\n";



open (SUPPREADS, ">", "$outprefix.supporting_reads.txt");
print SUPPREADS "chr\tstart\tend\tSTRlen\treadinfo\n";

# This used to be done per chromosome.  Changed by Joe.
print "Reading input data..\n";
open (F, "<", $inputfile) or die;
my $nread = 0;
my $indexPos = 0;
my @lenArray = ();
while (<F>) {
    ++$nread;
    if ($nread % 10000 == 0) {
        print "$nread lines processed\n";
    }

	chomp;
    my @element = split(" "); 
    my $chrom = $element[2];
    $chrom =~ s/chr//g;  # Get rid of 'chr' if it's there

    # Skip if not one of our requested chromosomes
    next if not exists $chrhash{$chrom};

	my $flag_direction = ""; # searching direction
	my $startRepeat = $element[3] + $element[10];
    my $endRepeat = $element[3] + $element[11]; 
    my $this_flank1 = substr($element[13], $element[10]-$flank_bp-1, $flank_bp);
    my $this_flank2 = substr($element[13], $element[11], $flank_bp);
	while (1) { # search all repeats for an overlapping record
		if ($endRepeat >= $repeatStart[$indexPos] and
            $startRepeat <= $repeatEnd[$indexPos]) {
			if ($element[9] eq $repeatType[$indexPos]) {
                # joe: fix the case where $element[9]-$flank_bp-1 is
                # negative, which causes the substring to be taken from
                # the end of the read
				if ($repeatFlank1[$indexPos] eq $this_flank1 and
				    $repeatFlank2[$indexPos] eq $this_flank2 and
                    $element[10] - $flank_bp - 1 >= 0) {
                    my $replen = scalar($element[11] - $element[10] + 1);
					if (exists $lenArray[$indexPos]) {
                        $lenArray[$indexPos] .= ",$replen";
                    } else {
                        $lenArray[$indexPos] .= "$replen";
                    }
                    print SUPPREADS "$chrom\t$repeatStart[$indexPos]\t$repeatEnd[$indexPos]\t$replen\t@element\n";
				}
			}
			last; # Don't allow a read to match multiple repeat records
		} elsif ($startRepeat > $repeatEnd[$indexPos]) { 
			# Forward searching (Genome repeat << Read repeat)
			if ($flag_direction eq "backward" or $indexPos >= $repeatno) {
                last;
            } else {
				$flag_direction = "forward";
                ++$indexPos;
            } 
		} elsif ($endRepeat < $repeatStart[$indexPos]) {
			# Backward searching (Read repeat >> Genome repeat)
			if ($flag_direction eq "forward" or $indexPos <= 0) {
                last;
            } else {
				$flag_direction = "backward";
                --$indexPos;
            }
		} else {
            print "ERROR: no repeat record found for read:\n";
			die("$chrom:$element[0]-$element[2]: $startRepeat-$endRepeat");
		}
	}
}
close F;
print "done.\n";

print "Writing results to $outprefix.str_summary.txt.. ";
open (OUTPUT, ">$outprefix.str_summary.txt");
print OUTPUT "index\tchr\tstart\tend\trepArray\n";
for (my $i = 0; $i < $repeatno; ++$i) {
    if (exists $lenArray[$i]) {
	    print OUTPUT "$i\t$repeatChrom[$i]\t$repeatStart[$i]\t$repeatEnd[$i]\t$lenArray[$i]\n";
    }
}
close OUTPUT;
print "done.\n";
