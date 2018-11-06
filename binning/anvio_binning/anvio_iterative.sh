# Daan Speth 2018
#
# A workflow for metagenome binning using anvio5 that relies on iterative assembly/binning
# rationale: 
# anvi'o interactive module can't handle too many (roughly >10k) contigs 
# taking a subset, bin contigs, repeat, allows manual binning with minimal hands on effort
#
# the workflow will:
# 1) (optional) map reads on pre-existing bins using bbmap to get unmapped reads
# 2) take a subset of 10 million (unmapped) reads (gives a number of contigs suitable for anvio binning; between 5k and 10k)
# 3) assemble using SPAdes
# 4) map 10M reads on SPAdes assembly, and select reads >95% identity over 80% of read, using BamM
# 5) Select contigs over 1000bp and 10x coverage for futher analysis
# 6) map all available datasets on selected contigs, and select reads >95% identity over 80% of read, using BamM
# 7) use anvi'o with ncbi COGs and centrifuge to build an annotated contig database
# 8) profile the contigs using all mapped datasets and merge profiles
#
# runs for 12-24 hrs on our server, and the result can be vizualized/binned using anvi-interactive:
# module load anaconda3 
# source activate anvio5
# anvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db
#
# bins can then be written using anvi-summarize (-C is collection_name)
# anvi-summarize -p SAMPLES-MERGED/PROFILE.db -c contigs.db -o dominant -C dominant
#
# written to work on a server with module system to manage $PATH variable 
# allows mixing of python2 (bamm) and python3 saftware (anvio)
#
# requires all data in a single datadir, I use softlinks
# assumes coupled paired end files (R1.fastq/R2.fastq)
# calls one custom perl script (also in this directory)


# set our defaults for the options
EXBIN_PATH=""
ANVI_COG_PATH="/data1/projects/dspeth/anvio/cog_collection/"
CENTRI_PATH="/data1/projects/dspeth/anvio/centifuge_db/"

# parse the options
while getopts 'e:c:a:' opt ; do
  case $opt in
    e) EXBIN_PATH=$OPTARG ;;
    c) ANVI_COG_PATH=$OPTARG ;;
    a) CENTRI_PATH=$OPTARG ;;
  esac
done
# skip over the processed options
shift $((OPTIND-1)) 

if [ $# -lt 3 ]
then
echo "Usage: $0 [options] iteration_number path_to_data_dir read_fastq_basename"
echo "Options: -e path_to_external_bins"
exit 1
fi

ITERATION="$1"	# number 
DATA="$2"		# absolute path to data-dir 
BASE="$3"		# fastq dataset basename (up until R1/R2)	

mkdir $BASE
cd $BASE
mkdir it_"$ITERATION"

if [ "$ITERATION" == "1" ] && [ -z "$EXBIN_PATH" ]
then
# get reads
cd it_"$ITERATION"
head -n 40000000 "$DATA"/"$BASE"_R1.fastq  > "$BASE"_R1_10Msubset.fastq 
head -n 40000000 "$DATA"/"$BASE"_R2.fastq  > "$BASE"_R2_10Msubset.fastq

# run spades
module load anaconda2
module load spades
spades.py --meta -o spades_10M -1 "$BASE"_R1_10Msubset.fastq -2 "$BASE"_R2_10Msubset.fastq -t 16
fi

if [ "$ITERATION" != "1" ] || [ -n "$EXBIN_PATH" ]
then
# get reads not assigned to previous bin
mkdir it_"$ITERATION"/unbinned 
cat it_*/anvio_out/dominant/bin_by_bin/Bin*/Bin*-contigs.fa > it_"$ITERATION"/unbinned/prev_it_bins_contigs.fa
cat "$EXBIN_PATH"/*fa > it_"$ITERATION"/unbinned/external_bins_contigs.fa
cat it_"$ITERATION"/unbinned/prev_it_bins_contigs.fa it_"$ITERATION"/unbinned/external_bins_contigs.fa > it_"$ITERATION"/unbinned/all_bins_co
ntigs.fa
cd it_"$ITERATION"
module load bbmap
bbmap.sh ref=unbinned/all_bins_contigs.fa in="$DATA"/"$BASE"_R1.fastq in2="$DATA"/"$BASE"_R2.fastq out=unbinned/"$BASE"_mapped.sam outu=unbin
ned/"$BASE"_unmapped.fastq threads=16 minid=0.90 idfilter=0.95 nodisk 
rm unbinned/"$BASE"_mapped.sam

# get reads - run spades
head -n 80000000 unbinned/"$BASE"_unmapped.fastq > "$BASE"_unmapped_10Msubset.fastq
module load spades
module load anaconda2
spades.py --meta -o spades_10M --12 "$BASE"_unmapped_10Msubset.fastq -t 16
fi

# map reads back
mkdir contig_selection
cp spades_10M/contigs.fasta contig_selection/contigs.fa
cd contig_selection
module load anaconda2 

if [ "$ITERATION" == "1" ] && [ -z "$EXBIN_PATH" ]
then
bamm make -d contigs.fa -c ../"$BASE"_R1_10Msubset.fastq ../"$BASE"_R2_10Msubset.fastq -t 16 -o raw
fi

if [ "$ITERATION" != "1" ] || [ -n "$EXBIN_PATH" ]
then
bamm make -d contigs.fa -i ../"$BASE"_unmapped_10Msubset.fastq -t 16 -o raw
fi


# Select contigs over 1000bp with coverage over 10
bamm filter -b raw/*bam -o id95_len80/ --percentage_id 0.95 --percentage_aln 0.80
bamm parse -c coverage.tsv -m opmean -b id95_len80/*bam
awk '$2 >= 1000 { print $0 }' coverage.tsv | awk '$3 >= 10 { print $0 }' > selected_contigs
blast_based_read_lookup_new.pl selected_contigs contigs.fa selected_contigs.fa

# prepare to map all reads to your selected contigs
cd ..
mkdir mapping
cd mapping
cp ../contig_selection/selected_contigs.fa .
module unload anaconda2 
module load anaconda3
source activate anvio5
anvi-script-reformat-fasta -o contigs_fixed.fa --simplify-names -r simpified_names selected_contigs.fa
perl -p -i -e s/\>/\>"$ITERATION"_/g contigs_fixed.fa
source deactivate
module unload anaconda3

# map all reads on the selected contigs, and filter bam file
module load anaconda2 
for i in "$DATA"/*R1.fastq ; do bamm make -d contigs_fixed.fa -c $i "${i:0:${#i}-8}"R2.fastq -t 16 -o raw; done
for i in raw/*bam ; do bamm filter -b $i -o id95_len80/ --percentage_id 0.95 --percentage_aln 0.80 ; done
module unload anaconda2

# anvio prep
cd ..
mkdir anvio_out
cp mapping/contigs_fixed.fa anvio_out/
cd anvio_out
module load anaconda3
source activate anvio5

# anvio contig db prep
anvi-gen-contigs-database -f contigs_fixed.fa -o contigs.db -n "$ITERATION"_"$BASE"
anvi-run-hmms -c contigs.db -T 16
anvi-run-ncbi-cogs -c contigs.db --cog-data-dir $ANVI_COG_PATH -T 16 
anvi-get-sequences-for-gene-calls -c contigs.db -o gene-calls.fa
centrifuge -f -x "$CENTRI_PATH"/p+h+v gene-calls.fa -S centrifuge_hits.tsv
anvi-import-taxonomy-for-genes -c contigs.db -i centrifuge_report.tsv centrifuge_hits.tsv -p centrifuge

# anvio profile
cd ../mapping/id95_len80/
for i in *bam; do anvi-profile -T 16 -M 1000 -S data_"${i:14:${#BASE}}" -o ../../anvio_out/"${i:14:${#BASE}}"_profile/ -c ../../anvio_out/con
tigs.db -i $i ; done
cd ../../anvio_out/
anvi-merge *profile/PROFILE.db -o SAMPLES-MERGED -c contigs.db --skip-concoct-binning -S sample_"$BASE"_"$ITERATION"
