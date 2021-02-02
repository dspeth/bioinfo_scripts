#! /usr/bin/perl
use strict;
use warnings;


# a quick and dirty script modified from a previously existing script to turn a protein alignment into SR4 numerical encoding
# If phylogeny software doesn't accept multi-state alignments, SR4_nuc_recode.pl can be used to "trick" the phylogeny software into accepting recoded alignment.
#
# input: alignment protein fasta file
#
# output: recoded alignment file
#
# Daan Speth, 2020

# requirements

if (scalar @ARGV != 2){
	die "\nRecodes a protein alignment to 4-state SR4 recoding in nucleotide form\nCan be used when phylogeny software supports multi-state alignments\n
  if multi-state alignments are not supported, SR4_nuc_recode.pl can be used\n\n
  requires two arguments:\n1) protein alignment fasta file\n2) user-defined outfile\n\nAGNPST -> 0\nCHWY -> 1\nDEKQR -> 2\nFILMV -> 3";
} 


# files

my $fasta = $ARGV[0];
my $out = $ARGV[1];

# read fasta into hash

my %sequences;
my $header = "placeholder";
open FASTA, $fasta or die "can not open fasta";
while (my $line = <FASTA>){
	chomp $line;
	my $first_char = substr($line,0,1);
	if ($first_char eq ">"){
		open OUT, (">> $out") or die "can not create file";
		if ( length $sequences{$header} ){
			$sequences{$header} =~ s/A/0/g; # SR4 group 1 AGNPST
                        $sequences{$header} =~ s/G/0/g;
			$sequences{$header} =~ s/N/0/g;
			$sequences{$header} =~ s/P/0/g;
			$sequences{$header} =~ s/S/0/g;
			$sequences{$header} =~ s/T/0/g;
			$sequences{$header} =~ s/C/1/g; # SR4 group 2 CHWY
                        $sequences{$header} =~ s/H/1/g;
			$sequences{$header} =~ s/W/1/g;
			$sequences{$header} =~ s/Y/1/g;
			$sequences{$header} =~ s/D/2/g; # SR4 group 3 DEKQR
			$sequences{$header} =~ s/E/2/g;
			$sequences{$header} =~ s/K/2/g;
			$sequences{$header} =~ s/Q/2/g;
			$sequences{$header} =~ s/R/2/g;
			$sequences{$header} =~ s/F/3/g; # SR4 group 4 FILMV 
			$sequences{$header} =~ s/I/3/g;
	                $sequences{$header} =~ s/L/3/g;
	                $sequences{$header} =~ s/M/3/g;
	                $sequences{$header} =~ s/V/3/g;
			print OUT "$sequences{$header}\n";
		}
                $header = substr($line, 1);
		print OUT ">$header\n";
		close OUT;
	}
	else {
		$sequences{$header} .= $line;
	}
	if (eof(FASTA)){
                        $sequences{$header} =~ s/A/0/g; # SR4 group 1 AGNPST
                        $sequences{$header} =~ s/G/0/g;
                        $sequences{$header} =~ s/N/0/g;
                        $sequences{$header} =~ s/P/0/g;
                        $sequences{$header} =~ s/S/0/g;
                        $sequences{$header} =~ s/T/0/g;
                        $sequences{$header} =~ s/C/1/g; # SR4 group 2 CHWY
                        $sequences{$header} =~ s/H/1/g;
                        $sequences{$header} =~ s/W/1/g;
                        $sequences{$header} =~ s/Y/1/g;
                        $sequences{$header} =~ s/D/2/g; # SR4 group 3 DEKQR
                        $sequences{$header} =~ s/E/2/g;
                        $sequences{$header} =~ s/K/2/g;
                        $sequences{$header} =~ s/Q/2/g;
                        $sequences{$header} =~ s/R/2/g;
                        $sequences{$header} =~ s/F/3/g; # SR4 group 4 FILMV 
                        $sequences{$header} =~ s/I/3/g;
                        $sequences{$header} =~ s/L/3/g;
                        $sequences{$header} =~ s/M/3/g;
                        $sequences{$header} =~ s/V/3/g;
		open OUT, ">> $out";
		print OUT "$sequences{$header}\n";
		close OUT;
	}

}



close FASTA;
