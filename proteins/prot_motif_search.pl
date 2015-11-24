#!/usr/bin/perl
use strict;
use warnings;

# a script to find a defined motif in a protein sequence.
# wildcards in the motif should be represented by a . 
# requires:


my $protein_fasta = $ARGV[0];
my $motif = uc $ARGV[1];
my $out_fasta = $ARGV[2]; 


if (@ARGV != 3){
	die "requires 3 arguments:\n\n1) protein fasta file\n2) user-defined motif\n3) user-defined outfile\n\n";
}


#read fasta into hash
my %prot;
my $prot_id;

open PROT, $protein_fasta or die "no protein fasta file provided";
while (my $line = <PROT>){
	chomp $line;
	my $fc = substr($line, 0, 1);
	if ($fc eq ">"){
		$prot_id = substr($line,1);
	}
	else{
		$prot{$prot_id} .= uc $line;
	}
}
close PROT;

# print out seqs that match input motif

foreach $prot_id (keys %prot){
	if ($prot{$prot_id} =~ /$motif/){
		close OUT;
    open OUT, ">> $out_fasta";
    print OUT ">$prot_id\n$prot{$prot_id}\n";
	}
}

