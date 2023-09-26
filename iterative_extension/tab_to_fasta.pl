#! /usr/bin/perl
use strict;
use warnings;

if (scalar @ARGV != 3){
	die "\nconvert a tab delimited file to a fasta file\n\nrequires 3 arguments\nusage: tab_file fasta_file column_containing_sequence\n\n";
}

my $contigs = $ARGV[0];
my $out = $ARGV[1];
my $seq_column = $ARGV[2];

my @stats;
my $count = 0;
 
open FILE, $contigs or die "can not find your file";
READ:	while (my $line = <FILE>){
		$count++ ;
		if ($count == 1){
			next READ;
		}
		else{
			chomp $line;
			my $fc = substr($line, 0, 1);
			@stats = split("\t", $line);
			my $header = $stats[0]; # change according to colum with header
			my $seq = $stats[$seq_column]; # change according to column with sequence
			open OUT,">> $out" or die "can not find file";
			print OUT "\>$header\n$seq\n";
			close OUT;
		}
	}
close FILE;
