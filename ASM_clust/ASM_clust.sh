#!/bin/bash -e

# set our defaults for the options
SUBSET=1000
PERP=1000
MAX_ITER=5000
ALIGN=mmseqs2
THREADS=1
FASTAREF=""

# parse the options
while getopts 's:p:m:a:t:f:' opt ; do
	case $opt in
		s) SUBSET=$OPTARG ;;
		p) PERP=$OPTARG ;;
		m) MAX_ITER=$OPTARG ;;
		a) ALIGN=$OPTARG ;;
		t) THREADS=$OPTARG ;;
		f) FASTAREF=$OPTARG ;;
	esac
done

# skip over the processed options
shift $((OPTIND-1))

if [[ "$ALIGN" != "mmseqs2" && "$ALIGN" != "diamond" && "$ALIGN" != "blast" ]]; then
	echo "allowed options for aligners (-a): mmseqs2, diamond, blastp"
	exit 1
fi

# check for mandatory positional parameters
if [ $# -lt 1 ]; then
	echo "Usage: $0 [options] fasta_file"
	echo "Options:"
	echo "-s INTEGER (subset of sequences for matrix, default 1000)"
	echo "-p INTEGER (t-SNE perplexity value, default 1000)"
	echo "-m INTEGER (max iterations of t-SNE, default 5000)"
	echo "-a mmseqs2/diamond/blast (aligner, one of three options, default mmseqs2)"
	echo "-t INTEGER (threads for aligner)"
	echo "-f FILE (fasta file of sequences to use as references)"
	exit 1
fi

hasCommand() {
	  command -v "$1" >/dev/null 2>&1 || { echo "Please make sure that $1 is in \$PATH."; exit 1; }
}

notExists "$1" && echo "$1 not found!" && exit 1;

# test whether required software is in path
hasCommand bhtsne.py
hasCommand blast_based_read_lookup_new.pl
hasCommand merge_tab_files_molyb2.pl

if [[ "$ALIGN" == "mmseqs2" ]]; then
	hasCommand mmseqs
elif [[ "$ALIGN" == "diamond" ]]; then
	hasCommand diamond
elif [[ "$ALIGN" == "blast" ]]; then
	hasCommand makeblastdb
	hasCommand blastp
fi

INPUT="$1"

#### check if first line of file starts with >
TEST=$(head -n 1 "$INPUT")
if [[ ${TEST:0:1} != ">" ]] ; then
	echo "input file does not apppear to be a fastA file"
       	exit 1
fi
if [ ! -z $FASTAREF ] ; then
	REFTEST=$(head -n 1 $FASTAREF)
		if [[ ${REFTEST:0:1} != ">" ]] ; then
        	echo "reference file does not apppear to be a fastA file"
        	exit 1
	fi
fi

#### check input file extension
if [ ${INPUT: -4} == ".faa" ]; then
	BASE=${INPUT:0:${#INPUT}-4}
elif [ ${INPUT: -3} == ".fa" ]; then
	BASE=${INPUT:0:${#INPUT}-3}
elif [ ${INPUT: -6} == ".fasta" ]; then
	BASE=${INPUT:0:${#INPUT}-6}
else
	echo "input file extension has to be .faa, .fa, or .fasta"
	exit 1
fi

##### pick subset

mkdir "$BASE"_"$ALIGN"
grep ">" "$INPUT" | awk '{print $1}' > "$BASE"_ids
perl -p -i -e 's/\>//g' "$BASE"_ids

if [ ! -z $FASTAREF ] ; then
	SUBSET=$(grep -c ">" "$FASTAREF")
	cp "$FASTAREF" "$BASE"_"$SUBSET"_ids.faa
	grep ">" "$FASTAREF" | awk '{print $1}' > "$BASE"_"$SUBSET"_ids
	perl -p -i -e 's/\>//g' "$BASE"_"$SUBSET"_ids
else
	shuf -n $SUBSET "$BASE"_ids > "$BASE"_"$SUBSET"_ids
	blast_based_read_lookup_new.pl "$BASE"_"$SUBSET"_ids $INPUT "$BASE"_"$SUBSET"_ids.faa
fi

mv "$BASE"_"$SUBSET"_ids.faa "$BASE"_"$SUBSET"_ids "$BASE"_ids "$BASE"_"$ALIGN"
cd "$BASE"_"$ALIGN"

##### Alignment section
# if aligning with mmseqs2
if [[ "$ALIGN" == "mmseqs2" ]]; then
	mmseqs easy-search ../"$INPUT" "$BASE"_"$SUBSET"_ids.faa "$BASE"_align tmp/ --threads $THREADS --max-seqs $SUBSET --mask 0
	rm -r tmp
fi

# if aligning with diamond
if [[ "$ALIGN" == "diamond" ]]; then
	diamond makedb --in "$BASE"_"$SUBSET"_ids.faa -d "$BASE"_"$SUBSET" -p $THREADS
	diamond blastp -q ../"$INPUT" -d "$BASE"_"$SUBSET" -o "$BASE"_align -p "$THREADS" -k "$SUBSET" --sensitive --masking no
fi

# if aligning with blast
if [[ "$ALIGN" == "blast" ]]; then
	makeblastdb -in "$BASE"_"$SUBSET"_ids.faa -dbtype prot
	blastp -query ../"$INPUT" -db "$BASE"_"$SUBSET"_ids.faa -num_threads $THREADS -evalue 0.001 -outfmt 6 -max_target_seqs $SUBSET -out "$BASE"_align
fi

# parse blast alignment file
while read id ; do grep -F "$id" "$BASE"_align > temp_"$id" ; done < "$BASE"_"$SUBSET"_ids
for i in temp_* ; do awk -v filt="${i:5}" 'BEGIN{FS=OFS="\t"} $2==filt {print $0}' $i > hits_"${i:5}" ; done
rm temp_*

for i in hits_* ; do merge_tab_files_molyb2.pl "$BASE"_ids "$i" score_"${i:5}" ; done
echo "seqID" >> "$BASE"_matrix
cat "$BASE"_ids >> "$BASE"_matrix
for i in score_* ; do cut -f 2 "$i" | paste "$BASE"_matrix - > "$BASE"_matrix_temp ; mv "$BASE"_matrix_temp "$BASE"_matrix ; done
cut -f 2- "$BASE"_matrix > "$BASE"_matrix_wo_id

##### Run TSNE
bhtsne.py -v -p $PERP -m $MAX_ITER -i "$BASE"_matrix_wo_id -o "$BASE"_matrix_wo_id_p"$PERP"_m"$MAX_ITER"
cp "$BASE"_matrix_wo_id_p"$PERP"_m"$MAX_ITER" ../"$BASE"_"$ALIGN"_tsne_p"$PERP"_m"$MAX_ITER"
