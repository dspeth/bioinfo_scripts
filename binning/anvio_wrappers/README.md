## wrapper shell scripts used in metagenome analyses with anvi'o

These two scripts can be used to facilitate manual binning or bin inspection with anvi'o. 
I originally wrote these when snakemake workflows weren't available for anvi'o, so their utility has likely been superseded.
That said, they have proven useful to me.


### requirements
These scripts were always written for specific use on the server I used myself, and will not work out of the box. 
They require an installation of conda, and environments called "assembly", "read_mapping", and "anvio_7".
Scripts use spades, coverm, samtools, bbmap, anvio, and a personal utility script "blast_based_read_lookup_new.pl" (included).


### anvio_map_reads_to_bin.sh
Shell script to map reads to a bin, and create anvio profiles. Allows for quick evaluation of bin quality. 
I use this for both my own MAGs and previously published bins I have an interest in.


### anvio_metagenome_subset.sh
Shell script to take a subset of 10 million reads of a metagenome, assemble with spades, map reads back, and create anvi'o profiles. 
I've found that using 10 million reads in assembly almost always results in a number of contigs with the range of the sweet spot for anvi'o visualization, and allows selective binning of highly abundant community members.
Can be run iteratively to take a 10M read subset after the initial binning. This usually doesn't result in many extra bins, but can give a few>
After running this script on the datasets of interest and manual binning, I typically perform a coassembly of all unbinned reads, bin with metabat, and then refine bins using anvio_map_reads_to_bin.sh


