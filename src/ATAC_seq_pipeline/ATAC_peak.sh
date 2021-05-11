#!/bin/bash
#Script to call peaks from BAMs from ATAC-seq.
#Run on individual sample; tools needed: samtools, bedtools, MACS2.
#3

#Effective genome size for build
GSIZE=2913022398
#Genebody bed file - from ChipPeakAnno R package
TSS=/storage1/fs1/ha/Active/maherlab/saha.d/projects/Data/GRCh38.TSS.bed
#Path to obtain bam file 
STORAGE=/storage1/fs1/ha/Active/maherlab/hdang/mCRPC-omics/atac-seq/alignments
#Blacklisted regions
BLACKLIST=/storage1/fs1/ha/Active/maherlab/saha.d/projects/Data/hg38-blacklist.v2.sorted.bed
#Sample ID
ID=$1
#Working dir for sample
WD=/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/mCRPC/"$ID"

cd "$WD"


#PEAK CALLING

#Call peaks with MACS in PE mode - for fragment length
macs2 callpeak -f BAMPE -g hs -n "$ID" -t "$WD"/"$ID".final.bam --outdir peaks 2> macs2.log


#Remove blacklist regions and non-chromosomal peaks
# bedtools intersect -v -a "$WD"/peaks/"$ID"_peaks.narrowPeak -b "$BLACKLIST"  | grep -P 'chr[0-9XY]+(?!_)' > "$WD"/"$ID"_filtered.peaks.narrowPeak


# #Calculate Frip scores
# TOTAL_READS=$(samtools view -c "$WD"/"$ID".final.bam)
# READS_IN_PEAKS=$(bedtools sort -i "$WD"/"$ID"_filtered.peaks.narrowPeak | bedtools merge -i stdin | bedtools intersect -u -abam "$WD"/"$ID".final.bam -b stdin -ubam -nonamecheck| samtools view -c )
# echo "$READS_IN_PEAKS/$TOTAL_READS" | bc -l > "$WD"/"$ID".frip.txt

