#!/usr/bin/perl
use strict;
use warnings;

# A script to extract reads that have BLAST/DIAMOND hits against a custom database
# assumes sequences in fastQ files are all single line

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


# read tab file containing BLAST/DIAMOND result against seq set of interest (DB) 
# write query ID and score into a hash for fast retrieval
my @hits;
my %score;
my $query;
my $prev_query = "muck_and_mock";

open BLAST, $blast_tab or die "no tabular blast file provided";
HIT: 	while (my $line = <BLAST>){
		chomp $line;
		@hits = split('\t', $line);
		$query = shift(@hits);
		if ($query eq $prev_query){
			next HIT;
		}
		else{
			$score{$query} = pop(@hits);
			$prev_query = $query;
		}
	}
close BLAST;	

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

#read fasta/fastq and print seq that match blast hits
my $read_id;
my @temp;
my $temp_id;
my $read_seq;
my $line_count = 0;
my $matches = 0;

open READ_FILE, $read_file or die "no read file provided";
SEQ: while (my $line = <READ_FILE>){
		chomp $line;
		$line_count++;
		if ($file_type eq "fastA"){
			my $fc = substr($line, 0, 1);
			if ($fc eq ">"){
				if ($matches == 1){ 
					close OUT;
					open OUT, ">> $out_fasta";
					print OUT "$read_seq\n";
					$matches = 0;					#reset $matches
					undef($read_seq);				#reset $read_seq
				}
				@temp = split(" ", $line);				# split header on spaces only take the first part, in case you're not blasting reads but sequences with longer headers
				$temp_id = $temp[0];
				$read_id = substr($temp_id, 1);
				if (exists $score{$read_id}){ 
					$matches = 1;
					close OUT;
					open OUT, ">> $out_fasta";
#					print OUT ">$read_id\n";			# prints just the ID (until the first space)	
					print OUT "$line\n"; 				# prints the entire header
				}
				else {
					next SEQ;
				}
			}
			else{
				if ($matches == 1){
					$read_seq .= $line;
					if (eof(READ_FILE)){
						close OUT;
						open OUT, ">> $out_fasta";
						print OUT "$read_seq\n";
						close OUT;
					}
				}
				else {
					next SEQ;
				}
			}
		}
		elsif ($file_type eq "fastQ"){
			if ($line_count == 1){		
				@temp = split(" ", $line);				# split header on spaces only take the first part, in case you're not blasting reads but sequences with longer headers
				$temp_id = $temp[0];
				$read_id = substr($temp_id, 1);
				if (exists $score{$read_id}){ 
					$matches = 1;
					close OUT;
					open OUT, ">> $out_fasta";
					print OUT ">$read_id\n";
				}
			}
			elsif ($line_count == 2){
				if ($matches == 1){ 
					close OUT;
					open OUT, ">> $out_fasta";
					print OUT "$line\n";
					$matches = 0;
				}
			$line_count = $line_count - 4;
			}
			else {
				next SEQ;
			}
		}
	}	
close READ_FILE;
