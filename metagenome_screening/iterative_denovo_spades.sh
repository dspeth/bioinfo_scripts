#!/bin/bash


#### a script to iteratively extend denovo assemblies using the SPAdes assembler and the BBmap (formerly bowtie2) aligner


#### it will then:
#### 1) assemble those into the seed sequence
#### 2) extract consensus
#### 3) map all previously extracted  reads from the metagenome using stringent parameters (0.98 identity)
#### 4) re-extract the mapped reads
#### 5) repeat the process for 15 cycles

#### parameter check
if [ $# -lt 2 ]
then
echo "Usage: $0 trusted_contigs reads"
exit 1
fi

CONTIGS="$1"
READS="$2"

### add required programs to PATH (only if using module system)
module load bbmap
module load spades

# 1) mapping
mkdir start_dir
mv $CONTIGS start_dir
cd start_dir/
bbmap.sh ref=$CONTIGS
bbmap.sh in=../$READS out=temp.sam idfilter=0.98 outm=reads_subset.fastq outputunmapped=f interleaved=true threads=20

rm temp.sam
SUBSET=reads_subset.fastq
echo $SUBSET

# organize
cp $CONTIGS $SUBSET ../
cd ../



################################################ iterative mapping ####################################################################

for i in {0..10..1}

do

WORKDIR=it_$i
mkdir $WORKDIR
mv $CONTIGS $SUBSET $WORKDIR

cd $WORKDIR

# assembly
SPADES_DIR=spades_$i

spades.py -o $SPADES_DIR --12 $SUBSET --careful --trusted-contigs $CONTIGS --tmp-dir /data1/tmp -k 21,33,55,77,99,127 -t 20

# organize
BBMAP_DIR=bbmap_$i
mkdir $BBMAP_DIR
mv $SPADES_DIR/contigs.fasta $BBMAP_DIR
cd $BBMAP_DIR
CONTIGS=contigs.fasta

#align all reads
bbmap.sh ref=$CONTIGS

bbmap.sh in=../../$READS out=temp.sam idfilter=0.98 outm=reads_subset.fastq outputunmapped=f interleaved=true threads=20

rm temp.sam
SUBSET=reads_subset.fastq

# organize
cp $CONTIGS $SUBSET ../../
cd ../../

done
