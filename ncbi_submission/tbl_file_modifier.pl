#! /usr/bin/perl

use strict;
use warnings;

# 10/10/2014 - Daan Speth - Radboud University Nijmegen

# a script to modify the .tbl file (prokka output) for submission to genbank
# IMPORTANT after running this script there will always be more manual corrections that have to be made to the .tbl file!

# this script is part of a workflow that I'm trying to work with, to facilitate WGS submission.
# 1) 	assemble (& bin) genome of interest

# 2) 	annotate with Prokka (Torsten Seemann - http://www.vicbioinformatics.com/software.prokka.shtml)
#		this script assumes you have disabled the "-c" flag in the prodigal line within prokka, to annotate partial genes

# 3)	If the genome is done using iontorrent sequencing, try to fix frameshift errors and reannotate. see github.com/dspeth/bioinfo_scripts/iontorrent errors

# 4) 	modify contig header of the prokka .fsa contig file as appropriate: [gcode=] [organism=] [strain=] [isolation-source=] 

# 5) 	create a template file at http://www.ncbi.nlm.nih.gov/WebSub/template.cgi

# 6) 	put template.sbt, .fsa contig file & .tbl annotation table in a single folder 

# 7) 	run tbl2asn (ncbi: http://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/) with recommended command: tbl2asn -t template.sbt -p . -M n -Z discrep

# 8) 	check error types in the errorsummary file  (http://www.ncbi.nlm.nih.gov/genbank/examples.wgs)

# 9) 	run this script and manually correct errors

# 10) 	check 'discrep' file for issues and if appropriate and general enough add to script

# requires a .tbl file and an .val file from the tbl2asn output

my $tbl = $ARGV[0];
my $val = $ARGV[1];
my $out	= $ARGV[2];

if (scalar @ARGV != 3){
	die "\nA script to modify ncbi annotation tables prior to submission\nuses .val file from tbl2asn output to correct .tbl file\n\nrequires 3 arguments\nusage: *.tbl *.val user-defined-outfile\n\n";
}	
# read .val file and extract the genes which have bad start and stop codons
# write into two seperate hashes

my %stop;
my %start;
my $cds_stop;
my $cds_start;
my $contig_length;

open ERRORS, $val or die "no error file provided";
while (my $line = <ERRORS>){
		chomp $line;
		if ($line =~ /NoStop/){
			($cds_stop) = ($line =~ /\-\> \[gnl\|Prokka\|(.+)\]/);		 # extracts gene identifier, make sure you annotated w Prokka in compliant mode
			($contig_length) = ($line =~ /len= (\d+)/); 
			$stop{$cds_stop} = $contig_length;
		}
		elsif ($line =~ /StartCodon/){
			($cds_start) = ($line =~ /\-\> \[gnl\|Prokka\|(.+)\]/); 	# extracts gene identifier, make sure you annotated w Prokka in compliant mode
			($contig_length) = ($line =~ /len= (\d+)/); 
			$start{$cds_start} = $contig_length;
		}
}
close ERRORS;

# correcting the annotation table

my @storage;
my $count = 0;

open TABLE, $tbl or die "no annotation table provided";
# some general corrections first
READ: while (my $line = <TABLE>){
		chomp $line;
		$line =~ s/\(\)//; 												# remove double brackets w/o content
		$line =~ s/\.$//; 												# remove period at the end of line
		$line =~ s/(gene\t)([A-Z]{3})/$1\L$2\E/;						# change case of gene names to NCBI appropriate
		$line =~ s/(gene\t)([A-Z])/$1\L$2\E/;							# change case of gene names to NCBI appropriate
		$line =~ s/proteins/protein/;									# getting rid of some plurals
		$line =~ s/subunits/subunit/;
		$line =~ s/\tunknown protein/\thypothetical protein/;			# saving NCBI the trouble
		$line =~ s/\:RefSeq\:/\:GenBank\:/;								# NCBI submission can't deal with the combo "AA sequence" & "RefSeq" in one 'inference' line			
		if ($line =~ /\tscore\t/){										# gets rid of the number of CRISPR repeats, since I can't find out how to annotate those
			next READ;
		}
		$line =~ s/product\t[Ss]imilar to/product\t/;					# some reference genomes (especially Kuenenia stuttgartiensis) contain a lot of 'similar to' annotations. tbl2asn changes all these to 'hypothetical protein', resulting in fatal errors
		$line =~ s/product\t[Ss]trongly similar to/product\t/;
		$line =~ s/ gene //;											# gets rid of the word gene from product names
		$line =~ s/\;/ \-/;												# 
# add more here as issues pile up
		
# now onto the fixing of start and stop codons 
		if ($line =~ /^\d/){											# acts on the lines defining the location of the feature
			open OUT, (">> $out") or die "outfile not specified";
			print OUT join("\n",@storage);
			print OUT "\n";
			close OUT;
			@storage = ();
			push(@storage, $line);
		}
		elsif ($line =~ /^>/){
			$count ++; 													# to make sure the file doesn't start with a blank line
			open OUT, (">> $out") or die "outfile not specified";
			print OUT join("\n",@storage);
			if ($count > 1){ 
				print OUT "\n";  
			}
			print OUT "$line";
			close OUT;
			@storage = ();
		}	
		elsif ($line =~ /locus_tag/){
			my $test = substr($line,-12); 								# extracts gene identifier, modify if identifyer is longer
			if (exists $start{$test}){									# partial at 5 prime end								
				my @split_line = split("\t",$storage[0]);	
				my $five_prime = $split_line[0];
				if ($storage[0] =~ /^\d\t/){							# codeblock to put in 5 prime partials at contig start 
					$storage[0] =~ s/^(\d)\t/\<1\t/;
					if (substr($storage[0], -4) !~ /gene/){ 
						push(@storage,"\t\t\tcodon_start\t$five_prime");
					}
				}
				else {													# codeblock to put in 5 prime partials at contig end 
					my $offset = ($start{$test} - $five_prime) + 1;			 
					$storage[0] =~ s/^\d+\t/\<$start{$test}\t/;
					if (substr($storage[0], -4) !~ /gene/){ 
						push(@storage,"\t\t\tcodon_start\t$offset");
					}						
				}
			}
			elsif(exists $stop{$test}){									# partial at 3 prime end
				my @split_line = split("\t",$storage[0]);	
				my $five_prime = $split_line[0];
				my $three_prime = $split_line[1];
				if ($storage[0] =~ /\t\d\t/){							# codeblock to put in 3 prime partials at contig start 
					$storage[0] =~ s/\t(\d)\t/\t\>1\t/;
				}
				else {													# codeblock to put in 3 prime partials at contig end 				 
					$storage[0] =~ s/\t\d+\t/\t\>$stop{$test}\t/;							
				}	
			}
			push(@storage, $line);
		}
		else { 
			push(@storage, $line);
		}
	}
close TABLE;

# don't forget to print out that last pesky gene feature
open OUT, (">> $out") or die "outfile not specified";
print OUT join("\n",@storage);
close OUT;

