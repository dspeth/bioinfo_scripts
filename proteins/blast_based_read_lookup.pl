#!/usr/bin/perl
use strict;
use warnings;

# A script to extract reads that have BLAST/DIAMOND hits against a custom database
# The script will read the entire read file into a hash, so might run into memory issues on enormous datasets


# requires:
# 1) tabular blast output, with the query in the first column (standard)
# 2) sequencing read fasta file
# 3) name of outfile (user defined)

my $blast_tab = $ARGV[0];
my $read_file = $ARGV[1];
my $out_fasta = $ARGV[2]; 

if (@ARGV != 3){
	die "requires 3 arguments:\n\n1) Tabular blast/diamond ouput file\n2) file containing sequencing reads used as BLAST/DIAMOND queries\n3) user-defined outfile\n\n";
}


# in case you put this script in a pipeline which has to deal both fasta and fastq files, determine filetype
my $file_type;

open READ_FILE, $read_file or die "no read file provided";
my $test = <READ_FILE>;
my $fc = substr($test, 0, 1);
if ($fc eq ">"){
	$file_type = "fastA";
}
elsif ($fc eq "@"){
	$file_type = "fastQ";
}
else{
	die "\nfile type of $read_file not recognized\ncheck file\n\n"
}

#read fasta/fastq into hash
my %reads;
my $read_id;
my @temp;
my $temp_id;
my $line_count = 0;

open READ_FILE, $read_file or die "no read file provided";
SEQ: while (my $line = <READ_FILE>){
		chomp $line;
		$line_count++;
		if ($file_type eq "fastA"){
			my $fc = substr($line, 0, 1);
			if ($fc eq ">"){
				@temp = split(" ", $line);				# split header on spaces only take the first part, in case you're not blasting reads but sequences with longer headers
				$temp_id = $temp[0];
				$read_id = substr($temp_id, 1);
			}
			else{
				$reads{$read_id} .= $line;
			}
		}
		elsif ($file_type eq "fastQ"){
			if ($line_count == 1){
				@temp = split(" ", $line);				# split header on spaces only take the first part, in case you're not blasting reads but sequences with longer headers
				$temp_id = $temp[0];
				$read_id = substr($temp_id, 1);
			}
			elsif ($line_count == 2){
				$reads{$read_id} = $line;
			$line_count = $line_count - 3;
			}
			else {
				next SEQ;
			}
		}
	}	
close READ_FILE;


# runs through blast file and prints sequences from hash
my @blast;
my $query_read;
my $prev_line = "start";

open BLAST, $blast_tab or die "no tabular blast file provided";
TAB:	while (my $line = <BLAST>){
			chomp $line;
			if ($line eq $prev_line){ 								# get rid of multiple hits in case you forgot to put only best hit
				next TAB;
			}
			@blast = split('\t', $line);
			$query_read = $blast[0];
			if (exists($reads{$query_read})){
				close OUT;
				open OUT, ">> $out_fasta";
				print OUT ">$query_read\n$reads{$query_read}\n";
			} 
		}
close BLAST;	
