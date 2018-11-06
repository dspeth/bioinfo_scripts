#! /usr/bin/perl
use strict;
use warnings;

if (scalar @ARGV != 2){
	die "\nconvert a tab delimited file to a fasta file\n\nrequires 2 arguments\nusage: tab_file fasta_file\n\n";
}

my $contigs = $ARGV[0];
my $out = $ARGV[1];

my @stats;
open FILE, $contigs or die "can not find your file";
while (my $line = <FILE>){
	chomp $line;
	my $fc = substr($line, 0, 1);
	@stats = split("\t", $line);
	my $header = $stats[0]; # change according to colum with header
	my $seq = $stats[4]; # change according to colum with sequence
	open OUT,">> $out" or die "can not find file";
	print OUT ">$header\n$seq\n";
	close OUT;
}
close FILE;
