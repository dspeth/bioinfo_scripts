#! /bin/bash
eval "$(conda shell.bash hook)"
conda deactivate

# shell script to extend contigs obtained from a metagenome assembly, using the source reads
# useful for connecting low coverage contigs, works by reassembling a subset of reads mapped to the contigs
# PacBio mapping code added to the script by Bernhard Tschitschko

# set defaults for the options, in this case names of conda envs containing the required software
ASSEMBLY="assembly"
MAP="read_mapping"

while getopts 'a:m:' opt ; do
    case $opt in
        a) ASSEMBLY=$OPTARG ;;
        m) MAP=$OPTARG ;;
    esac
done

# skip over the processed options
shift $((OPTIND-1))

#### parameter check
if [ $# -lt 7 ]
then
echo "Usage: $0 seed_contigs_path pacbio_path fw_reads_path rv_reads_path project_name no_of_iterations threads"
echo ""
echo "requires 7 positional arguments:"
echo "contigs:		path to fasta file with seed contigs"
echo "pacbio reads:		path to the reverse reads"
echo "fw reads:		path to the forward reads"
echo "rv reads:		path to the reverse reads"
echo "project name:	prefix for files"
echo "num iterations: 	number of iterations to run the loop"
echo "num threads: 	number of threads to use where appropriate"
echo ""
echo "optional arguments:"
echo "-a name of conda env with read mappers, default read_mapping"
echo "-m name of conda env with spades, default assembly"
exit 1
fi


SEED=$1
PACBIO=$2
FW=$3
RV=$4
NAME=$5
END=$6
THREADS=$7

conda deactivate
mkdir it_1
cd it_1

conda activate $MAP

## previous version of the script used bbmap, current one uses minimap. Old command:
# bbmap.sh ref=$SEED in=$FW in2=$RV nodisk outm="$NAME"_mapped.sam threads=$THREADS minid=0.9 > bbmap_stdout 2>&1

# map PacBio and filter pacbio HiFi reads
minimap2 -ax map-hifi --sam-hit-only -t $THREADS $SEED $PACBIO > pacbio_vs_${NAME}_aln.sam
samtools view -bS -@ 15 pacbio_vs_${NAME}_aln.sam > pacbio_vs_${NAME}_aln.bam
coverm filter -b pacbio_vs_${NAME}_aln.bam -o pacbio_vs_${NAME}_filtered.bam -t 40 --min-read-percent-identity 98
samtools fasta -0 pacbio_vs_${NAME}_filtered.fasta -o pacbio_vs_${NAME}_filtered.fasta -1 pacbio_vs_${NAME}_filtered.fasta -2 pacbio_vs_${NAME}_filtered.fasta pacbio_vs_${NAME}_filtered.bam

# map Illumina short reads
minimap2 -t $THREADS -a -x sr --score-N 2 --sam-hit-only --no-pair $SEED $FW $RV > ${NAME}_illuminamap.sam 2> minimap_stderr
samtools view -bShu ${NAME}_illuminamap.sam | samtools sort -m 97G -@ 4 - -o ${NAME}_illuminamap_sorted.bam
samtools index ${NAME}_illuminamap_sorted.bam

# filter mapped reads 98% id, 50% mapped
coverm filter --bam-files ${NAME}_illuminamap_sorted.bam -o ${NAME}_illuminamap_sorted_filtered.bam --min-read-percent-identity 98 --min-read-aligned-percent 50
samtools sort -o ${NAME}_illuminamap_sorted_filtered_sort.bam --threads $THREADS ${NAME}_illuminamap_sorted_filtered.bam
samtools fastq -1 ${NAME}_illuminamap_98_80_R1.fastq -2 ${NAME}_illuminamap_98_80_R2.fastq ${NAME}_illuminamap_sorted_filtered_sort.bam

# complete broken pairs 
seqkit seq -n ${NAME}_illuminamap_98_80_R1.fastq ${NAME}_illuminamap_98_80_R2.fastq | sort | uniq > ${NAME}_illuminamap_98_80_ids
seqkit grep -f ${NAME}_illuminamap_98_80_ids $FW > ${NAME}_illuminamap_98_80_pair1.fastq
seqkit grep -f ${NAME}_illuminamap_98_80_ids $RV > ${NAME}_illuminamap_98_80_pair2.fastq

conda deactivate

# assemble using spades
conda activate $ASSEMBLY

# hybrid assembly
spades.py --isolate -t $THREADS -k 21,33,55,77,99,111 -s pacbio_vs_${NAME}_filtered.fasta -1 ${NAME}_illuminamap_98_80_pair1.fastq -2 ${NAME}_illuminamap_98_80_pair2.fastq --trusted-contigs $SEED -o spades_${NAME} > spades_stdout 2>&1

conda deactivate

# old way of filtering for length, could/should be replaced by seqkit command
cd spades_${NAME}
fasta_to_gc_cov_length_tab.pl scaffolds.fasta scaf.tab
awk '$3 >= 2000 {print $0}' scaf.tab > scaf_filtered.tab
tab_to_fasta.pl scaf_filtered.tab scaf_filtered.fa 3
cd ../../

# loops over the rest of the iterations
START=2
for ITERATION in $(eval echo "{$START..$END}") ; do
	mkdir it_"$ITERATION"
	cp it_$(("$ITERATION"-1))/spades_"$NAME"/scaf_filtered.fa it_"$ITERATION"
	cd it_"$ITERATION"
	conda activate $MAP

	minimap2 -ax map-hifi --sam-hit-only -t $THREADS scaf_filtered.fa $PACBIO > pacbio_vs_${NAME}_aln.sam
	samtools view -bS -@ 15 pacbio_vs_${NAME}_aln.sam > pacbio_vs_${NAME}_aln.bam
	coverm filter -b pacbio_vs_${NAME}_aln.bam -o pacbio_vs_${NAME}_filtered.bam -t 40 --min-read-percent-identity 98
	samtools fasta -0 pacbio_vs_${NAME}_filtered.fasta -o pacbio_vs_${NAME}_filtered.fasta -1 pacbio_vs_${NAME}_filtered.fasta -2 pacbio_vs_${NAME}_filtered.fasta pacbio_vs_${NAME}_filtered.bam


 	minimap2 -t $THREADS -a -x sr --score-N 2 --sam-hit-only --no-pair scaf_filtered.fa $FW $RV > ${NAME}_illuminamap.sam 2> minimap_stderr
	samtools view -bShu ${NAME}_illuminamap.sam | samtools sort -m 97G -@ 3 - -o ${NAME}_illuminamap_sorted.bam
	samtools index ${NAME}_illuminamap_sorted.bam


	coverm filter --bam-files ${NAME}_illuminamap_sorted.bam -o ${NAME}_illuminamap_sorted_filtered.bam --min-read-percent-identity 98 --min-read-aligned-percent 80
	samtools sort -o ${NAME}_illuminamap_sorted_filtered_sort.bam --threads $THREADS ${NAME}_illuminamap_sorted_filtered.bam
	samtools fastq -1 ${NAME}_illuminamap_98_80_R1.fastq -2 ${NAME}_illuminamap_98_80_R2.fastq ${NAME}_illuminamap_sorted_filtered_sort.bam

  	seqkit seq -n ${NAME}_illuminamap_98_80_R1.fastq ${NAME}_illuminamap_98_80_R2.fastq | sort | uniq > ${NAME}_illuminamap_98_80_ids
  	seqkit grep -f ${NAME}_illuminamap_98_80_ids $FW > ${NAME}_illuminamap_98_80_pair1.fastq
  	seqkit grep -f ${NAME}_illuminamap_98_80_ids $RV > ${NAME}_illuminamap_98_80_pair2.fastq

	conda deactivate

	conda activate $ASSEMBLY
	spades.py --isolate -t $THREADS -k 21,33,55,77,99,111 -s pacbio_vs_${NAME}_filtered.fasta -1 ${NAME}_illuminamap_98_80_pair1.fastq -2 ${NAME}_illuminamap_98_80_pair2.fastq --trusted-contigs scaf_filtered.fa -o spades_${NAME} > spades_stdout 2>&1
  	conda deactivate

	# old way of filtering for length, could/should be replaced by seqkit command
  	cd spades_${NAME}
	fasta_to_gc_cov_length_tab.pl scaffolds.fasta scaf.tab
	awk '$3 >= 2000 {print $0}' scaf.tab > scaf_filtered.tab
	tab_to_fasta.pl scaf_filtered.tab scaf_filtered.fa 3
	cd ../../
done
