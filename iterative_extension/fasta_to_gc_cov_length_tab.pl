#! /usr/bin/perl
use strict;
use warnings;


# a script to turn contig fasta files exported from the clc genomics workbench in to a tab delimited form usable in R
# since the contig header includes the average coverage after assembly it becomes one of the columns
# If you are parsing a fasta file from a different source, just change the printed header accordingly (line 48)
#
# 
# input: fasta file
#
# output: tab delimited list with these columns:
# Contig_name	Sequencing_depth	GC_content	Contig_length	Sequence
# (WARNING: each line contains the full sequence, can be a bit of a pain when checking file with 'head')	
#
# Daan Speth, 2014 

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
# SPAdes assembler output
#		$line =~ s/_cov_/\t/;
# CLC assembler output
#		$line =~ s/ Average coverage\: /\t/; # this specifcally splits the ID output of the CLC genomics workbench into ID and the value for depth  
		$line =~ s/,//; #remove commas from contigs with depth over 1000 to avoid pesky R issues
		$header = substr($line, 1);
		}
	else {
		$sequences{$header} .= $line;
	}
}
close FASTA;

# print hash id and gc content
my $total;
my $gc;
my $gc_content;
open OUT, (">> $out") or die "can not create file";
# 	print OUT "Contig\tSequencing depth\tGC content\tContig length\tSequence\n"; #adjust line according to headers. Suggestions below
# 	print OUT "Contig\tSequencing depth\tGC content\tContig length\n";
	print OUT "Contig\tGC content\tContig length\tSequence\n";
foreach $header(keys(%sequences)) {
	$total = length($sequences{$header});
	$gc = ($sequences{$header} =~ tr/GCgc//);
	$gc_content = $gc / $total;
# either of the lines below can be used to print the info in tabular format
# the longer one also prints the sequences. 
# comment out the one you don't need
	print OUT "$header\t$gc_content\t$total\t$sequences{$header}\n";
#	print OUT "$header\t$gc_content\t$total\n";
}
close OUT;
