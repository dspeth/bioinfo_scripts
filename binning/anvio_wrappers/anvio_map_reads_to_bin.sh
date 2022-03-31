#!/bin/bash

# a script to map data on an anvio bin and display for further analysis 

# set our defaults for the options
INTERLEAVED=0
THREADS=1

while getopts 'it:' opt ; do
	case $opt in
		i) INTERLEAVED=1 ;;
		t) THREADS=$OPTARG ;;
	esac
done

# skip over the processed options
shift $((OPTIND-1))

if [ $# -lt 5 ]
then
echo "Usage: $0 contigs data_dir bin_name outdir fastq_name_length"
echo ""
echo "contigs: the absolute path to contigs fasta file"
echo "data dir: the absolute path to the data dir containing fastq files"
echo "bin name: the chosen name for the bin/MAG"
echo "out_dir: name of a directory that will be created for the output files"
echo "fastq_name_length: the length of the fastq filename up to _val_1/_val_2"
echo ""
echo "optional arguments:"
echo "specify number of threads using \"-t NUM\" (default 1)"
echo "when fastq files are interleaved (and not paired) use \"-i\""
echo ""
echo "the script assumes the name format of paired-end fastq files is *val_1.fq.gz and *val_2.fq.gz"
echo "this is the format of trim_galore trimmed files, with a specified prefix"
echo "if the filenames differ from this format, please edit to comply"
echo ""
exit 1
fi

CONTIGS="$1"	# absolute path to contigs
DATA="$2"		# absolute path to data-dir
NAME="$3" 		# bin/genome name
OUT="$4"		# output directory
LENGTH="$5"		# length of fastq basename, for distinguising anvio profiles

# anvio prep
mkdir $OUT
cd $OUT
mkdir "$NAME"_anvio_out
cd "$NAME"_anvio_out
eval "$(conda shell.bash hook)"
conda activate anvio-7.1

# anvio contig db prep
anvi-script-reformat-fasta -o "$NAME"_reformat.fasta --simplify-names -r contig_name_report $CONTIGS
anvi-gen-contigs-database -f "$NAME"_reformat.fasta -o contigs.db -n $NAME
anvi-run-hmms -c contigs.db -T $THREADS
anvi-run-ncbi-cogs -c contigs.db -T $THREADS
anvi-run-pfams -c contigs.db -T $THREADS
anvi-run-scg-taxonomy -c contigs.db -T $THREADS
anvi-run-kegg-kofams -c contigs.db -T $THREADS

# prepare to map all reads to your contigs
cd ..
mkdir "$NAME"_mapping
cd "$NAME"_mapping

# map all reads on the contigs, and filter bam file

conda deactivate
conda activate read_mapping
# line for paired files
if [[ $INTERLEAVED == 0 ]]; then
	for i in "$DATA"/*1.fq.gz ; do coverm make -r ../"$NAME"_anvio_out/"$NAME"_reformat.fasta --coupled $i "${i:0:${#i}-7}"2.fq.gz -t $THREADS -o raw; done
fi
# line for interleaved files
if [[ $INTERLEAVED == 1 ]]; then
	for i in "$DATA"/*fastq ; do coverm make -r ../"$NAME"_anvio_out/"$NAME"_reformat.fasta --interleaved $i -t $THREADS -o raw; done
fi

mkdir id95_len80
cd raw
for i in *bam ; do samtools sort -@ $THREADS -o sort_"$i" $i ; done
for i in sort*bam ; do coverm filter -b $i -o ../id95_len80/"$i" --min-read-percent-identity 95 --min-read-aligned-percent 80 ; done
conda deactivate

# anvio profile
cd ../id95_len80/
conda activate anvio-7.1
for i in *bam ; do anvi-init-bam -o anvio_"$i" -T $THREADS $i ; done
for i in anvio_*bam ; do anvi-profile -T $THREADS -M 1000 -S sample_"${i:${#NAME}+27:$LENGTH}" -o ../../"$NAME"_anvio_out/sample_"${i:${#NAME}+27:$LENGTH}"_profile/ -c ../../"$NAME"_anvio_out/contigs.db -i $i ; done
cd ../../"$NAME"_anvio_out/
anvi-merge *profile/PROFILE.db -o SAMPLES-MERGED -c contigs.db -S "$NAME"
