#!/usr/bin/perl
use strict;
use warnings;

# a script to correct errors in assemblies based on iontorrent reads
# Iontorrent data suffers from strand specific errors, which can persist in assembly, resulting in frameshifts. 
# If you care about correcting these, this script might be of use. I doubt it, however, since it is tailored to my own needs. 
# Feel free to modify though

# the script uses the outcome of a BLAST and a variant caller to insert nucleotides in positions probably causing frameshifts

# a suggested workflow (which is typically a pain):

# 1) annotate with prokka (torsten seemann) 
# 2) map reads, call variants (both clc) (settings probabilistic variant caller: frequency >30; count > 10; fwd/rev < 0.3; avg qual >15; length = 1 type = insertion
# 3) blast proteins to reference to find frameshifts (clc or commandline as long as the output is a tab delimited file, containing the queryID (in order of the genome and the hitID)
#	 check below if lines 63-69 are suitable
# 4) check the blast file
# 5) use script below to correct errors broken genes
# repeat annotation and manually inspect broken genes 

# requires: 
# contigs in fasta format
# test files (list of proteins with frameshifts, blast file)
# - blast file
# - tab separated file with the location of insertions containing columns: contig id, type (Insertion or Deletion), position, nucleotide, annotated feature disrupted
# also, this indel file should be ordered by indel position, from furthest from the beginning to nearest (n to 1)

# input format should be made less restrictive, but haven thought of a way to do so yet, suggestions welcome @ d.speth at science dot ru dot nl :)

# set input parameters
my $draft_contigs = $ARGV[0];
my $blast_tab = $ARGV[1];
my $indel = $ARGV[2];
my $blast_out = $ARGV[3];
my $out = $ARGV[4];

#load fasta into hash
my %contigs;
my $header;

open DRAFT, $draft_contigs or die "no contig file provided";
while (my $line = <DRAFT>){
	chomp $line;
	my $fc = substr($line, 0, 1);
	if ($fc eq ">"){
		$header = substr($line,1);
	}
	else{
		$contigs{$header} .= $line;
	}
}
close DRAFT;

my %broken_orfs;
my $test;
my @blast;
my $ORF;
my $hit;
my $prev_ORF = 'startvalue';
my $prev_hit = 'startvalue';

open BLAST, $blast_tab or die "no tabular blast file provided";
READ: while (my $line = <BLAST>){
		chomp $line;
		@blast = split('\t', $line);

##### uncomment lines when using a CLC blast table where the hit is in the last column ######		
		$ORF = shift(@blast);
		$hit = pop(@blast);
#############################################################################################

#		$ORF = shift(@blast);
#		$hit = shift(@blast);
		if ($ORF eq $prev_ORF){ ### only take the best hit, in case you forget to put the max no of hits at 1...
			next READ;
		}
		if ($hit eq $prev_hit){
			$test = $prev_ORF;
			$broken_orfs{$test} = $prev_hit;
			$test = $ORF;
			$broken_orfs{$test} = $hit;
		}	
		$prev_ORF = $ORF;
		$prev_hit = $hit;	
}
close BLAST;

# sanity check
my $n_broken_orfs = (scalar keys(%broken_orfs)) / 2;
print "number of broken ORFS:\n$n_broken_orfs\n";

open (BLAST_OUT, ">> $blast_out") or die "cannot create file";
foreach my $test (keys %broken_orfs){
	print BLAST_OUT "$test\t$broken_orfs{$test}\n"; 
}


# read through indel table and insert nucleotides where appropriate
my $insert_count = 0;
my $del_count = 0;
my @indel_feats;
my @processed_ORFS;
print "\nORF with multiple corrected positions\ncheck these ORFS manually and consider rerunning script\n";

open INDEL, $indel or die "no indels found, probably deleted, insert please...";	
while (my $modifier = <INDEL>){
	chomp $modifier;
	@indel_feats = split("\t", $modifier);
	my $real_errors = pop(@indel_feats);
	if (exists($broken_orfs{$real_errors})){	
		if ($real_errors ~~ @processed_ORFS){
			print "$real_errors\n";
		}
		push(@processed_ORFS, $real_errors);
		my $id = shift(@indel_feats);
		my $position = shift(@indel_feats);
		my $type = shift(@indel_feats);
		my $nucleotide = shift(@indel_feats);
		if (exists($contigs{$id})){
			if ($type eq 'Insertion'){
				substr($contigs{$id}, $position, 0) = "$nucleotide";
				$insert_count ++;
			}
			elsif ($type eq 'Deletion'){
				$position = $position - 1;
				substr($contigs{$id}, $position, 1, "");
				$del_count ++;
			}
		}
	}
}
close INDEL;
	
# print modified contig hash in fasta format
open (OUTPUT, ">> $out") or die "cannot create file";
foreach my $id (keys %contigs){
	print OUTPUT ">$id\n$contigs{$id}\n"; 
}
print "...\n\ncontigs modified\n\n$insert_count nucleotides inserted\n$del_count nucleotides deleted\n";

			
