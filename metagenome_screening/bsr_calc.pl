#!/usr/bin/perl
use strict;
use warnings;

# requires:
# 1) tabular blast output against gene of interest, with the query in the first column (standard)
# 2) tabular blast output against outgroup, with the query in the first column (standard)
# 3) name of outfile (user defined)

my $blast_db = $ARGV[0];
my $blast_og = $ARGV[1];
my $bsr_tab = $ARGV[2]; 

if (@ARGV != 3){
	die "requires 3 arguments:\n\n1) Tabular blast/diamond ouput against gene of interest\n2) Tabular blast/diamond ouput against outgroup \n3) user-defined outfile to store tabular BSR \n\n";
}


# read tab file containing BLAST/DIAMOND result against seq set of interest (DB) 
# write query ID and score into a hash for fast retrieval
my @hits_db;
my $query_db;
my %score_db;
my %length_db;
my %id_db;
my %hit_start_db;
my %hit_id_db;
my $prev_q_db = "mock_n_muck";

open BLAST_DB, $blast_db or die "no tabular blast file provided";
DB:	while (my $line = <BLAST_DB>){
		chomp $line;
		@hits_db = split('\t', $line);
		$query_db = $hits_db[0];
		if ($query_db eq $prev_q_db){
			next DB;
		}
		else {
			$score_db{$query_db} = pop(@hits_db);
			$length_db{$query_db} = $hits_db[3];
			$id_db{$query_db} = $hits_db[2];
			$hit_start_db{$query_db} = $hits_db[8];
			$hit_id_db{$query_db} = $hits_db[1];
		}
		$prev_q_db = $query_db;
}
close BLAST_DB;

# read tab file containing BLAST/DIAMOND result against outgroup (OG) 
# outgroup should include seq set of interest as well
# write query ID and score into a hash for fast retrieval
my @hits_og;
my $query_og;
my %score_og;
my $prev_q_og = "muck_and_mock";

open BLAST_OG, $blast_og or die "no tabular blast file provided";
OG:	while (my $line = <BLAST_OG>){
		chomp $line;
		@hits_og = split('\t', $line);
		$query_og = $hits_og[0];
		if ($query_og eq $prev_q_og){
			next OG;
		}
		else{
			$score_og{$query_og} = pop(@hits_og);
		}
		$prev_q_og = $query_og;
}
close BLAST_OG;

# calculate bit score ratio (BSR) and print file

my $bsr;

open OUT, ">> $bsr_tab";
print OUT "query\thit\tdb_score\tog_score\tbsr\taln_length\tpct_id\thit_start\n";
close OUT;

foreach my $id (keys(%score_og)){ 
	if (exists($score_db{$id})){
		$bsr = $score_db{$id} / $score_og{$id};
		close OUT;
		open OUT, ">> $bsr_tab";
		print OUT "$id\t$hit_id_db{$id}\t$score_db{$id}\t$score_og{$id}\t$bsr\t$length_db{$id}\t$id_db{$id}\t$hit_start_db{$id}\n";
	} 
}


# print file with DB matches that don't have a hit in the outgroup db

my $no_hit = $bsr_tab . "_no_outgroup_hits";

open NOHIT, ">> $no_hit";
print NOHIT "query\thit\tdb_score\taln_length\tpct_id\thit_start\n";
close NOHIT;

foreach my $id (keys(%score_db)){
	if (! exists($score_og{$id})){
		open NOHIT, ">> $no_hit";                                                                               # due to the differences in database size (very small fo
r gene of interest, large$
		print NOHIT "$id\t$hit_id_db{$id}\t$score_db{$id}\t$length_db{$id}\t$id_db{$id}\t$hit_start_db{$id}\n";	# these are logged for later inspection, as they may inc
lude divergent hits to gene of interest.
		close NOHIT;
	}
}

close NOHIT;
close OUT;	

