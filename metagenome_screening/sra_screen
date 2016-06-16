#! /bin/bash

#### sra_screen ####
#
# A script to systematically download a predefined set of metagenomes and search them for the presence of a marker gene of interest.
# uses wget, sra-toolkit, DIAMOND and a custom perl script (www.github.com/dspeth)
# requires making a DIAMOND database prior to running script
# 
# output for each dataset searched: diamond files and a read fasta file containing the reads of interest
#
# Daan Speth, 2016
#
####################

#### parameter check
if [ $# -lt 3 ]
then
echo "Usage: $0 SRA-accession-list DIAMOND-db processed-accession-list"
exit 1
fi

SRA="$1"
DMND_DB="$2"
PROCESSED="$3"

# load the required modules
module load sra_toolkit

# make directories for output
mkdir dmnd_out fasta_files

# read through file with SRA accessions and download/process each file in succession 
while read ACCESSION 
do 
	echo "processing $ACCESSION" >> log.txt
	
	# Check if Accession number has already been processes
	if grep -Fxq $ACCESSION $PROCESSED
	then
		echo "$ACCESSION already processed - Skipping" >> log.txt
		continue	
	fi

	# download dataset
	wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/"${ACCESSION:0:3}"/"${ACCESSION:0:6}"/"$ACCESSION"/"$ACCESSION".sra 

	# unpack dataset from .sra format into .fastq format & split reads
	fastq-dump --split-3 "$ACCESSION".sra
	rm "$ACCESSION".sra 
	
	# because of proccessing time, we only look at one of the split files in case of paired reads
	if [ -e "$ACCESSION"_1.fastq ]
	then
		rm "$ACCESSION".fastq "$ACCESSION"_2.fastq
		diamond blastx -q "$ACCESSION"_1.fastq -d $DMND_DB -k 1 -a hits_"$ACCESSION" --seg no -t /export/data1/tmp/
	# for single reads, data will be put i a single file
	else 	
		diamond blastx -q "$ACCESSION".fastq -d $DMND_DB -k 1 -a hits_"$ACCESSION" --seg no -t /export/data1/tmp/
	fi
	
	# convert DIAMOND outfile to tab delimited format
	diamond view -a hits_"$ACCESSION" -o tab_"$ACCESSION"
	
	# use tab delimited file for read lookup in the original file using custom perl script
	if [ -e "$ACCESSION"_1.fastq ]
	then
		blast_based_read_lookup.pl tab_"$ACCESSION" "$ACCESSION"_1.fastq seqs_"$ACCESSION".fasta
		rm "$ACCESSION"_1.fastq
	else 
		blast_based_read_lookup.pl tab_"$ACCESSION" "$ACCESSION".fastq seqs_"$ACCESSION".fasta
		rm "$ACCESSION".fastq
	fi
	
	echo "$ACCESSION" >> $PROCESSED
	
	mv hits_"$ACCESSION".daa tab_"$ACCESSION" dmnd_out/			
	mv seqs_"$ACCESSION".fasta fasta_files/

done < $SRA
