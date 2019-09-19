#! /usr/bin/perl
use strict;
use warnings;

# a script to merge two tab delimited alignment files by the values in their first column.
# specifically for use with the ASM_clust.sh shell script

# if used for other purposes, be sure to pay attention to the following

# script presumes the score (column of interest) is the last column
# script presumes all sequences are aligned to a single reference sequence
# script presumes the 2nd file is a subset of the first and fills in 0 for values which are not present in the 2nd file.

# Daan Speth 2019

if (scalar @ARGV != 3){
        die "\nmerge two tab delimited files by their first column\n\nrequires three arguments:\n1) 1st tab delimited file\n2) 2nd tab delimited file\n3) user-defined temp outfile (can overwrite infile)\n\n";
}


# files

my $tab_1 = $ARGV[0];
my $tab_2 = $ARGV[1];
my $out = $ARGV[2];

# read second table file into hash
my @addition;
my %score_2;
my $first_2;
my $ref_2;

open ADD, $tab_2 or die "can not open second file";
while (my $line = <ADD>){
	chomp $line;
	@addition = split("\t", $line);	
	$first_2 = shift @addition;
	$ref_2 = shift @addition;
	$score_2{$first_2} = pop @addition;
}
close ADD;

my @columns;
my $first_1;
my $count = 1;

open OUT, ("> $out") or die "can not create file";
open BASE, $tab_1 or die "can not create file";

while (my $line = <BASE>){
        chomp $line;
	@columns = split("\t", $line);
	$first_1 = shift @columns;
	if ($count == 1){
		print OUT "seqID\t$ref_2\n";
		$count ++;
	}
	if (exists($score_2{$first_1})){
		print OUT "$line\t$score_2{$first_1}\n";
        }
	else{
		print OUT "$line\t0\n";
	}
}

close BASE;
close OUT;
