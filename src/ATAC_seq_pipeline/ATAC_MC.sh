#!/bin/bash
#Script to assess significance of read overlap at TSS for hg38 using Monte Carlo simulation via permutation testing

#TSS bed file of protein coding tx to check enrichment from Ensdb R package
TSS=/storage1/fs1/ha/Active/maherlab/saha.d/projects/Data/TSS_ensdb.v86_hg38.bed 
#Dir to store files
STORAGE=/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/Files
#Chromosome size txt from UCSC
CHROM=/storage1/fs1/ha/Active/maherlab/saha.d/projects/Data/hg38.chrom.sizes
#SampleID
ID=$1
#PE Beds for each sample
BEDPE="$STORAGE"/Final_bedpe/"$ID".bedpe.bed
#Sorted bedfile for reads in sample
BED="$STORAGE"/"$ID"_filter_sort_reads.bed

bedtools sort -i "$BEDPE" > "$BED"

#Find number of intersections of PE bed with TSS region - Append unshuffled intersections to start of output
bedtools intersect -u -a "$BED" -b "$TSS" -nonamecheck -sorted | wc -l >>  Permutation_counts/"$ID".perm_counts.txt

#For 100 iterations, permute TSS locations and find intersections - append to output file
for i in `seq 1 100`;
do  
    bedtools shuffle -i "$TSS" -g "$CHROM" -chrom | bedtools sort -i stdin | bedtools intersect -u -a "$BED" -b stdin -nonamecheck -sorted | wc -l >>  Permutation_counts/"$ID".perm_counts.txt
done



