The whole genome FASTA used for making the bed file of regions containing the CTCF motif is not included because of its size. It can be found through UCSC, and decompressed with gzip using the following two commands:

wget http://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz

The generate_exact-match-motif_bed_CTCF.py script expects hg38.fa to be found in this (helper_scripts) directory.
