#!/usr/bin/env bash

# workflow for metagenome analysis using anvio
# relies on iterative assembly/binning
# requires all data in a single datadir
# assumes coupled paired end files (R1.fastq/R2.fastq)


# set our defaults for the options
EXBIN_PATH=""

# parse the options
while getopts 'e:' opt ; do
  case $opt in
    e) EXBIN_PATH=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1))

if [ $# -lt 4 ]
then
echo "Usage: $0 [options] iteration_number path_to_data_dir read_fastq_basename threads"
echo "Options: -e path_to_external_bins"
exit 1
fi

ITERATION="$1"	# integer
DATA="$2"		# absolute path to data-dir
BASE="$3"		# fastq dataset basename (up until R1/R2)
THREADS="$4"	# integer, number of threads used for spades/bamm

if [[ "$BASE" =~ [^a-zA-Z0-9_] ]]; then
echo "the basename of the fastq files should only contain alphanumeric characters or an underscore"
echo "please change basename of all the files to contain only alphanumeric characters and underscores"
exit 1
fi

mkdir $BASE
cd $BASE
mkdir it_"$ITERATION"
eval "$(conda shell.bash hook)"

if [ "$ITERATION" == "1" ] && [ -z "$EXBIN_PATH" ]
then
# get reads
cd it_"$ITERATION"
head -n 40000000 "$DATA"/"$BASE"_R1.fastq  > "$BASE"_R1_10Msubset.fastq 
head -n 40000000 "$DATA"/"$BASE"_R2.fastq  > "$BASE"_R2_10Msubset.fastq

# run spades
conda activate assembly
spades.py --meta -o spades_10M -1 "$BASE"_R1_10Msubset.fastq -2 "$BASE"_R2_10Msubset.fastq -t $THREADS
conda deactivate
fi

if [ "$ITERATION" != "1" ] || [ -n "$EXBIN_PATH" ]
then
# get reads not assigned to previous bin
mkdir it_"$ITERATION"/unbinned
cat it_*/anvio_out/dominant/bin_by_bin/Bin*/Bin*-contigs.fa > it_"$ITERATION"/unbinned/prev_it_bins_contigs.fa
cat "$EXBIN_PATH"/*fa > it_"$ITERATION"/unbinned/external_bins_contigs.fa
cat it_"$ITERATION"/unbinned/prev_it_bins_contigs.fa it_"$ITERATION"/unbinned/external_bins_contigs.fa > it_"$ITERATION"/unbinned/all_bins_contigs.fa
cd it_"$ITERATION"

conda activate read_mapping
bbmap.sh ref=unbinned/all_bins_contigs.fa in="$DATA"/"$BASE"_R1.fastq in2="$DATA"/"$BASE"_R2.fastq out=unbinned/"$BASE"_mapped.sam outu=unbinned/"$BASE"_unmapped.fastq threads=$THREADS minid=0.90 idfilter=0.95 nodisk
conda deactivate

rm unbinned/"$BASE"_mapped.sam

# get reads - run spades
head -n 80000000 unbinned/"$BASE"_unmapped.fastq > "$BASE"_unmapped_10Msubset.fastq
conda activate assembly
spades.py --meta -o spades_10M --12 "$BASE"_unmapped_10Msubset.fastq -t $THREADS
conda deactivate
fi

# map reads back
mkdir contig_selection
cp spades_10M/contigs.fasta contig_selection/contigs.fa
cd contig_selection

conda activate read_mapping
if [ "$ITERATION" == "1" ] && [ -z "$EXBIN_PATH" ]
then
coverm contig -r contigs.fa --coupled ../"$BASE"_R1_10Msubset.fastq ../"$BASE"_R2_10Msubset.fastq -o coverage.tab --min-read-percent-identity 95 --min-read-aligned-percent 80 -t $THREADS
fi

if [ "$ITERATION" != "1" ] || [ -n "$EXBIN_PATH" ]
then
coverm contig -r contigs.fa --interleaved ../"$BASE"_unmapped_10Msubset.fastq -o coverage.tab --min-read-percent-identity 95 --min-read-aligned-percent 80 -t $THREADS
fi
conda deactivate

# Select contigs over 1000bp with coverage over 10

fasta_to_gc_cov_length_tab.pl contigs.fa contigs_gc_len.tab
sort -k1 -o coverage.tab coverage.tab
sort -k1 -o contigs_gc_len.tab contigs_gc_len.tab
join -t $'\t' -1 1 -2 1 -o 1.1,1.2,1.3,2.2 contigs_gc_len.tab coverage.tab > contigs_gc_len_cov.tab
awk '$3 >= 1000 && $4 >= 10 {print $0}' contigs_gc_len_cov.tab > selected_contigs
blast_based_read_lookup_new.pl selected_contigs contigs.fa selected_contigs.fa

# prepare to map all reads to your selected contigs
cd ..
mkdir mapping
cd mapping
cp ../contig_selection/selected_contigs.fa .

conda activate anvio7
anvi-script-reformat-fasta -o contigs_fixed.fa --simplify-names -r simpified_names selected_contigs.fa
perl -p -i -e s/\>/\>"$ITERATION"_/g contigs_fixed.fa
conda deactivate

# map all reads on the selected contigs, and filter bam file
conda activate read_mapping
for i in "$DATA"/*R1.fastq ; do coverm make -r contigs_fixed.fa --coupled $i "${i:0:${#i}-8}"R2.fastq -t $THREADS -o raw; done
cd raw
for i in *bam ; do samtools sort -@ $THREADS -o sort_"$i" $i ; done
mkdir ../id95_len80/
for i in sort*bam ; do coverm filter -b $i -o ../id95_len80/"$i" --min-read-percent-identity 95 --min-read-aligned-percent 80 ; done
conda deactivate


# anvio prep
cd ../..
mkdir anvio_out
cp mapping/contigs_fixed.fa anvio_out/
cd anvio_out
conda activate anvio7

# anvio contig db prep
anvi-gen-contigs-database -f contigs_fixed.fa -o contigs.db -n "$ITERATION"_"$BASE"
anvi-run-hmms -c contigs.db -T $THREADS
anvi-run-ncbi-cogs -c contigs.db -T $THREADS
anvi-run-pfams -c contigs.db -T $THREADS
anvi-run-scg-taxonomy -c contigs.db -T $THREADS
anvi-run-kegg-kofams -c contigs.db -T $THREADS

# anvio profile
cd ../mapping/id95_len80/
for i in *bam ; do anvi-init-bam -o "${i:0:${#i}-4}"_init.bam $i ; done
for i in *init.bam; do anvi-profile -T $THREADS -M 1000 -S data_"${i:22:${#BASE}}" -o ../../anvio_out/"${i:22:${#BASE}}"_profile/ -c ../../anvio_out/contigs.db -i $i ; done
cd ../../anvio_out/
anvi-merge *profile/PROFILE.db -o SAMPLES-MERGED -c contigs.db -S sample_"$BASE"_"$ITERATION"
