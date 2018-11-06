#! /usr/bin/perl
use strict;
use warnings;

# a script to assign contigs to bins of the tetraESOM pipeline

# I have used initial ESOM binning to visually guide differential coverage binning
# first run ESOM as described here: https://github.com/tetramerfreqs/binning
# then bin the data and run this script.
# merge the bin info with other data in R and use the bin assignment as fill in scatterplot. Contigs not assgned to bins will have no fill


# needs 4 arguments:
# 1) .cls file containing the numbers of binned sequences coupled to their bins
# 2) .names files containing the numbers of sequences coupled to their headers
# 3) .fasta file of original contigs used as ESOM input
# 4) bin_list: name for a new file containing the contig headers and and the bin designation 
# 5) prefix for the bin files 

####### IMPORTANT: make sure fasta headers match the headers of the contigs after ESOM (names file)
####### (all perl incompatible characters are substituted by underscores or removed)


# input arguments
my $bin = $ARGV[0];
my $names = $ARGV[1];
my $contigs = $ARGV[2];
my $bin_list = $ARGV[3];
my $out = $ARGV[4];

if (scalar @ARGV != 5){
die "\nA script to correct indels in iontorrent files\n\nrequires 5 arguments\nusage: ESOM_cls_file ESOM_names_file contig_fasta bin_list bin_prefix\n\n";
}


# open bin file and write into hash
my %bins;
my $id;
my @pair;
open BIN, $bin or die "can not find bin file";
BINNING: while (my $line = <BIN>){
		chomp $line;
		my $first_char = substr($line, 0, 1);
		if ($first_char eq "%"){
			next BINNING;
		}
# the following is probably a retarded way to read two values into a hash, but I can't come up with a better solution right now... (and it works)
		else {
			@pair = split(/\t/, $line);
			$id = shift @pair;
			$bins{$id} = pop @pair;
		}
}
close BIN;

# open name file and write lines into array and write hash of name + class
my %classes;
my $header;
my $index;
my @names;
open NAMES, $names or die "what names file are you talking about???";
CLASSIFY: while (my $line = <NAMES>) { 
			chomp $line;
			my $first_char = substr($line, 0, 1);
			if ($first_char eq "%"){
				next CLASSIFY;
			}
# the following code will create an array of the line, take the last value (the contig identifier) as header of a new hash and couple the bin nr to the contig identifier			
			else {
				@names = split(/\t/, $line);
				$header = pop @names;
				$index = shift @names;
				if (exists($bins{$index})){
					$classes{$header} = $bins{$index};
				}
			}
}
close NAMES;
undef %bins; # frees up some memory

#prints the bin designation coupled to header for uses in plots later
open BIN_LIST, ">> $bin_list" or die "can't write list";
foreach $header (keys %classes){
	print BIN_LIST "$header\t$classes{$header}\n";
	}
close BIN;

# read through fasta and print
my $contig_id; 
open CONTIGS, $contigs or die "can not find contigs";
while (my $line = <CONTIGS>){
	chomp $line;
	$contig_id = substr($line,1);
	if (exists($classes{$contig_id})){
		close OUT;
		open OUT, ">> $out\_$classes{$contig_id}\.fa";
		print OUT ">bin_$classes{$contig_id}\_$contig_id\n";
	} 
	else {
		print OUT "$line\n";
	}
}
close OUT;

