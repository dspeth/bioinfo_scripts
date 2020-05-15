### ASM clust scripts have been moved to their own repo at https://github.com/dspeth/ASM_clust
### please use those versions, rather than the ones in this directory



ASM-clust

the three scripts in this directory comprise ASM-clust, an approach to de novo classify complex protein superfamilies
ASM-clust is intended to 

ASM-Clust is implemented in bash with two helper scripts in perl, and will take a protein fasta file as the sole input. 
Fasta files are processed with ASM_clust.sh, which then:
1) randomly selects a subset of n sequences (default 1000) 
2) aligns the entire dataset to the subset of n sequences
3) combines all scores into a matrix (inserting 0 for query-database pairs that did not produce an alignment)
4) reduces the matrix to 2 dimensions using t-SNE (Van der Maaten and Hinton 2008; Van der Maaten 2014)

External dependencies for ASM-clust are:
1) the python wrapper for the Barnes-Hut implementation of t-SNE available here :(https://github.com/lvdmaaten/bhtsne)
2) Alignment software. For flexible usage, ASM-Clust supports alignment using:
    1) MMSeqs2 (Steinegger and Söding 2017) (default aligner)
    2) DIAMOND (Buchfink, Xie, and Huson 2015)
    3) BLAST (Altschul et al. 1990)

Other user-defined options are:
- the number of sequences in the subset (default 1000)
- the main t-SNE parameter “perplexity” (default 1000) 
- t-SNE maximum iterations (default 5000) 
- the number of threads used by the alignment software (default 1). 

Although the clustering is generally similar with multiple randomly chosen subsets, 
the subset can be defined for reproducibility. 

The output of ASM_clust.sh can be visualized as a scatterplot where each dot represents a sequence, 
and clusters are readily apparent. This format allows additional annotation with sequence features, 
such as taxonomy, length, or composition. 
