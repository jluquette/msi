# TODO:
#  make flank_size a command line arg
#  use a better perl argument parsing library
#  fix the field indeces in the parsing at the end of the script.  (i.e.,
#  since I added the mapQ value to the line, everything needs to shift over
#  by one.)

chomp $ARGV[0];
$filename = $ARGV[0];
chomp $ARGV[1];
$realfile = $ARGV[1];

@chrarray = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
#@chrarray = qw(10 Y);

print "$filename ..\n";
open (output0, ">SPUTNIK_".$filename."_".$chrarray[0].".txt");
open (output1, ">SPUTNIK_".$filename."_".$chrarray[1].".txt");
open (output2, ">SPUTNIK_".$filename."_".$chrarray[2].".txt");
open (output3, ">SPUTNIK_".$filename."_".$chrarray[3].".txt");
open (output4, ">SPUTNIK_".$filename."_".$chrarray[4].".txt");
open (output5, ">SPUTNIK_".$filename."_".$chrarray[5].".txt");
open (output6, ">SPUTNIK_".$filename."_".$chrarray[6].".txt");
open (output7, ">SPUTNIK_".$filename."_".$chrarray[7].".txt");
open (output8, ">SPUTNIK_".$filename."_".$chrarray[8].".txt");
open (output9, ">SPUTNIK_".$filename."_".$chrarray[9].".txt");
open (output10, ">SPUTNIK_".$filename."_".$chrarray[10].".txt");
open (output11, ">SPUTNIK_".$filename."_".$chrarray[11].".txt");
open (output12, ">SPUTNIK_".$filename."_".$chrarray[12].".txt");
open (output13, ">SPUTNIK_".$filename."_".$chrarray[13].".txt");
open (output14, ">SPUTNIK_".$filename."_".$chrarray[14].".txt");
open (output15, ">SPUTNIK_".$filename."_".$chrarray[15].".txt");
open (output16, ">SPUTNIK_".$filename."_".$chrarray[16].".txt");
open (output17, ">SPUTNIK_".$filename."_".$chrarray[17].".txt");
open (output18, ">SPUTNIK_".$filename."_".$chrarray[18].".txt");
open (output19, ">SPUTNIK_".$filename."_".$chrarray[19].".txt");
open (output20, ">SPUTNIK_".$filename."_".$chrarray[20].".txt");
open (output21, ">SPUTNIK_".$filename."_".$chrarray[21].".txt");
open (outputX, ">SPUTNIK_".$filename."_".$chrarray[22].".txt");
open (outputY, ">SPUTNIK_".$filename."_".$chrarray[23].".txt");
#open (data, "SPUTNIKresult_exome_len5_".$filename.".txt") or die;
open (data, $realfile) or die;
while ($line = <data>) {
@element = split(" ", $line); 
$element[2] =~ s/chr//g; $element[2] = "chr".$element[2];
if ($element[2] eq "chr".$chrarray[0]) {
	print output0 $line;
} elsif ($element[2] eq "chr".$chrarray[1]) {
	print output1 $line;
} elsif ($element[2] eq "chr".$chrarray[2]) {
	print output2 $line;
} elsif ($element[2] eq "chr".$chrarray[3]) {
	print output3 $line;
} elsif ($element[2] eq "chr".$chrarray[4]) {
	print output4 $line;
} elsif ($element[2] eq "chr".$chrarray[5]) {
	print output5 $line;
} elsif ($element[2] eq "chr".$chrarray[6]) {
	print output6 $line;
} elsif ($element[2] eq "chr".$chrarray[7]) {
	print output7 $line;
} elsif ($element[2] eq "chr".$chrarray[8]) {
	print output8 $line;
} elsif ($element[2] eq "chr".$chrarray[9]) {
	print output9 $line;
} elsif ($element[2] eq "chr".$chrarray[10]) {
	print output10 $line;
} elsif ($element[2] eq "chr".$chrarray[11]) {
	print output11 $line;
} elsif ($element[2] eq "chr".$chrarray[12]) {
	print output12 $line;
} elsif ($element[2] eq "chr".$chrarray[13]) {
	print output13 $line;
} elsif ($element[2] eq "chr".$chrarray[14]) {
	print output14 $line;
} elsif ($element[2] eq "chr".$chrarray[15]) {
	print output15 $line;
} elsif ($element[2] eq "chr".$chrarray[16]) {
	print output16 $line;
} elsif ($element[2] eq "chr".$chrarray[17]) {
	print output17 $line;
} elsif ($element[2] eq "chr".$chrarray[18]) {
	print output18 $line;
} elsif ($element[2] eq "chr".$chrarray[19]) {
	print output19 $line;
} elsif ($element[2] eq "chr".$chrarray[20]) {
	print output20 $line;
} elsif ($element[2] eq "chr".$chrarray[21]) {
	print output21 $line;
} elsif ($element[2] eq "chr".$chrarray[22]) {
	print outputX $line;
} elsif ($element[2] eq "chr".$chrarray[23]) {
	print outputY $line;
} 
}
close data; close output0; close output1; close output2; close output3; close output4;
close output5; close output6; close output7; close output8; close output9;
close output10; close output11; close output12; close output13; close output14;
close output15; close output16; close output17; close output18; close output19;
close output20; close output21; close outputX; close outputY;


# Read chromosome fasta sequences
@chrarray = qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X Y);
#@chrarray = qw(10 Y);
$chrSeq = (); 
foreach $chr (@chrarray) {
	print "reading chromosome $chr .. \n";
	open (data, "/scratch/ljl11/neuron/msitools/hg19_chromosome_seq/chr".$chr.".fa") or die;
	$null = <data>;	while ($line = <data>) {chomp $line; $chrSeq{"chr".$chr} .= $line;}
	close data; 
}

open (output, ">SUMMARY_wgsCalls_".$filename.".txt");
print output "index\tchr\tstart\tend\trepArray\n";
open (SUPPREADS, ">SUPPORT_reads_".$filename.".txt");
print SUPPREADS "chr\tstart\tend\tSTRlen\treadinfo\n";

$flank_size = 10;

foreach $chr (@chrarray) {
	@repeatStart = (); @repeatEnd = (); @repeatSeq = (); @repeatFlank1 = (); @repeatFlank2 = (); @repeatType = (); $repeatno = 0;
	open (data, "/scratch/ljl11/neuron/msitools/WGRef_7892585MS_withGENEcategory_FINAL_withMSTYPE_hg19.txt") or die;
	while ($line = <data>) {
		chomp $line; @element = split("\t", $line); 
		if ($element[0] eq "chr".$chr) {
			$element[5] =~ s/\s//g;
			$repeatType[$repeatno] = $element[5]; $repeatStart[$repeatno] = $element[1]; 
			$repeatEnd[$repeatno] = $element[2];	$repeatSeq[$repeatno] = $element[4];
			$repeatFlank1[$repeatno] = uc(substr($chrSeq{"chr".$chr}, $element[1] - ($flank_size+1), $flank_size));
			$repeatFlank2[$repeatno] = uc(substr($chrSeq{"chr".$chr}, $element[2] , $flank_size));
			++$repeatno;
		}
	}
	close data;

	print "Chr $chr - $filename - WGS .. \n";	@foundRead = (); 
	open (data, "SPUTNIK_".$filename."_".$chr.".txt") or die;
	$indexPos = 0; @lenArray = (); $flag_direction = "forward";
	while (<data>) {
		chomp; @element = split(" "); 
		$startRepeat = $element[3] + $element[9]; $endRepeat = $element[3] + $element[10]; 
		$flag_direction = ""; # searching direction
		while (1) { # searching overlapping repeat from indexPos
			#print "$indexPos - $flag_direction for $startRepeat to $endRepeat VS $repeatStart[$indexPos] to $repeatEnd[$indexPos] "; $null = <STDIN>; 

			if ($endRepeat >= $repeatStart[$indexPos] and $startRepeat <= $repeatEnd[$indexPos]) {
				if ($element[8] eq $repeatType[$indexPos]) {
	                #print "flank1 : $repeatFlank1[$indexPos] and flank2 : $repeatFlank2[$indexPos] \n";

                #print substr($element[12],$element[9]-($flank_size+1),$flank_size)."::".substr($element[12],$element[10],$flank_size)."\n";
					#print "read: ".substr($element[12],$element[9]-($flank_size+1),$element[10]-$element[9]+4)."\n"; # $null = <STDIN>;
                    # joe: fix the case where $element[9]-3 is negative, which causes the substring to be taken from the end of the read
					if ($repeatFlank1[$indexPos] eq substr($element[12],$element[9]-($flank_size+1),$flank_size) and $element[9]-($flank_size+1) >= 0 and $repeatFlank2[$indexPos] eq substr($element[12],$element[10],$flank_size)) {
						$lenArray[$indexPos] .= scalar($element[10]-$element[9]+1).",";
                        print SUPPREADS "$chr\t$repeatStart[$indexPos]\t$repeatEnd[$indexPos]\t".scalar($element[10]-$element[9]+1)."\t@element\n";
					}
				}
				last; # terminate since a read is allowed for a 'single' match
			} elsif ($startRepeat > $repeatEnd[$indexPos]) { 
				# Forward searching (Genome repeat << Read repeat) ----------------------------------
				if ($flag_direction eq "backward" or $indexPos >= $repeatno) {last} else {
					$flag_direction = "forward"; ++$indexPos;} 
			} elsif ($endRepeat < $repeatStart[$indexPos]) {
					# Backward searching (Read repeat >> Genome repeat) ---------------------------------
				if ($flag_direction eq "forward" or $indexPos <= 0) {last} else {
					$flag_direction = "backward"; --$indexPos;}
			} else {
				print "Failed! READ: $element[0] - $element[2] : $startRepeat - $endRepeat"; exit;
			}
			#if ($indexPos < 0) {++$indexPos} elsif ($indexPos >= $repeatno) {last}
		}
	}
	close data;
	for ($i = 0; $i < $repeatno; ++$i) {
		print output $i."\t$chr\t$repeatStart[$i]\t$repeatEnd[$i]\t$lenArray[$i]\n";
	}
	#unlink "SPUTNIK_".$filename."_".$chr.".txt";
} # End of chromosome loop
close output;

