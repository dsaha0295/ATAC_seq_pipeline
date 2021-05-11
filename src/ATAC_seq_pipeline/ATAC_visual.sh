#!/bin/bash
#Script to visualize peaks on BAMs from ATAC-seq. 
#Run on individual sample; tools needed: deeptools
#2

#Effective genome size for build
GSIZE=2913022398
#TSS bed file - from ChipPeakAnno R package
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

#PLOTS

#Create bigwig files for visualization from coordinate sorted filtered BAMs (RPGC = reads/bin normalized to 1x coverage)
bamCoverage -b "$WD"/"$ID".final.bam -o "$WD"/"$ID".bw  --binSize 10 --normalizeUsing RPGC --effectiveGenomeSize "$GSIZE" 

#TSS enrichment plot
computeMatrix reference-point -S "$WD"/"$ID".bw -R "$TSS"  --referencePoint TSS  -a 2000 -b 2000 -out "$WD"/"$ID"_TSS.tab.gz 

plotProfile -m "$WD"/"$ID"_TSS.tab.gz -out "$WD"/"$ID".png --plotTitle "$ID TSS enrichment"

rm "$WD"/"$ID"_TSS.tab.gz 
