#! /usr/bin/perl
use strict;
use warnings;


# a script to turn a fasta file into a tab delimited file with columns "ID" and "Length"
#
# Daan Speth, 2016 

# requirements

if (scalar @ARGV != 2){
	die "\ncalculate contig stats and return tab delimited file, including sequence\n\nrequires two arguments:\n1) fasta file\n2) user-defined outfile\n\ncheck script for format specifics\n";
} 

# files

my $fasta = $ARGV[0];
my $out = $ARGV[1];

# read fasta into hash

my %sequences;
my $header;
open FASTA, $fasta or die "can not open fasta";
while (my $line = <FASTA>){
	chomp $line;
	$line =~ s/\r//;
	my $first_char = substr($line,0,1);
	if ($first_char eq ">"){
		$header = substr($line, 1);
		}
	else {
		$sequences{$header} .= $line;
	}
}
close FASTA;

# print hash id and length
my $seqlength;
open OUT, (">> $out") or die "can not create file";
	print OUT "Protein_id\tLength\n";
foreach $header(keys(%sequences)) {
	$seqlength = length($sequences{$header});
	print OUT "$header\t$seqlength\n";
}
close OUT;
