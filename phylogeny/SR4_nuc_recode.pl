#! /usr/bin/perl
use strict;
use warnings;


# a quick and dirty script modified from a previously existing script to turn a protein alignment into SR4 mucleotide encoding
# Usable to "trick" phylogeny software into accepting a recoded alignment as nucleotides.
# If phylogeny software accept multi state alignments (e.g. RAxML), use SR4_numerical_recode.pl
#
# input: alignment fasta file
#
# output: recoded alignment file
#
# Daan Speth, 2020

# requirements

if (scalar @ARGV != 2){
	die "\nRecodes a protein alignment to 4-state SR4 recoding in nucleotide form\nCan be used when phylogeny software doesn't support multi-state alignments\n
  Of multi-state alignments are supported (e.g. RAxML), use SR4_numerical_recode.pl\n\n
  requires two arguments:\n1) fasta file\n2) user-defined outfile\n\nAGNPST -> A\nCHWY -> C\nDEKQR -> G\nFILMV -> T";
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
			$sequences{$header} =~ s/G/A/g; # SR4 group 1 AGNPST
			$sequences{$header} =~ s/N/A/g;
			$sequences{$header} =~ s/P/A/g;
			$sequences{$header} =~ s/S/A/g;
			$sequences{$header} =~ s/T/A/g;
			$sequences{$header} =~ s/H/C/g; # SR4 group 2 CHWY
			$sequences{$header} =~ s/W/C/g;
			$sequences{$header} =~ s/Y/C/g;
			$sequences{$header} =~ s/D/G/g; # SR4 group 3 DEKQR
			$sequences{$header} =~ s/E/G/g;
			$sequences{$header} =~ s/K/G/g;
			$sequences{$header} =~ s/Q/G/g;
			$sequences{$header} =~ s/R/G/g;
			$sequences{$header} =~ s/F/T/g; # SR4 group 4 FILMV 
			$sequences{$header} =~ s/I/T/g;
	    $sequences{$header} =~ s/L/T/g;
	    $sequences{$header} =~ s/M/T/g;
	    $sequences{$header} =~ s/V/T/g;
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
                $sequences{$header} =~ s/G/A/g; # SR4 group 1 AGNPST
                $sequences{$header} =~ s/N/A/g;
                $sequences{$header} =~ s/P/A/g;
                $sequences{$header} =~ s/S/A/g;
                $sequences{$header} =~ s/T/A/g;
                $sequences{$header} =~ s/H/C/g; # SR4 group 2 CHWY
                $sequences{$header} =~ s/W/C/g;
                $sequences{$header} =~ s/Y/C/g;
                $sequences{$header} =~ s/D/G/g; # SR4 group 3 DEKQR
                $sequences{$header} =~ s/E/G/g;
                $sequences{$header} =~ s/K/G/g;
                $sequences{$header} =~ s/Q/G/g;
                $sequences{$header} =~ s/R/G/g;
                $sequences{$header} =~ s/F/T/g; # SR4 group 4 FILMV 
                $sequences{$header} =~ s/I/T/g;
                $sequences{$header} =~ s/L/T/g;
                $sequences{$header} =~ s/M/T/g;
                $sequences{$header} =~ s/V/T/g;
		open OUT, ">> $out";
		print OUT "$sequences{$header}\n";
		close OUT;
	}

}


close FASTA;
