#! usr/bin/perl

use strict;
use warnings;

my $fasta_file = $ARGV[0];
my $tab_file = $ARGV[1];

if (scalar @ARGV != 2){
	die "\nconvert a fasta file to a tab delimited files\n\nrequires 2 arguments\nusage: fasta_file tab_file\n\n";
}

# fasta to tab

#load fasta into hash
my %fasta;
my $header;

open DRAFT, $fasta_file or die "no contig file provided";
while (my $line = <DRAFT>){
	chomp $line;
	my $fc = substr($line, 0, 1);
	if ($fc eq ">"){
		$header = substr($line,1);
	}
	else{
		$fasta{$header} .= $line;
	}
}
close DRAFT;

# print modified contig hash in fasta format
open (OUTPUT, ">> $tab_file") or die "cannot create file";
foreach $header(keys %fasta){
	print OUTPUT ">$header\t$fasta{$header}\n"; 
}
close OUTPUT;
