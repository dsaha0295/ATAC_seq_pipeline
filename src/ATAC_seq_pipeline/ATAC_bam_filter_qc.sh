#!/bin/bash
#Script to filter and perform QC on BAMs from ATAC-seq - based on ENCODE and Harvard FAS informatics pipelines.
#Run on individual sample; tools needed: samtools, bedtools, Picard.
#1

#Path to obtain bam file 
STORAGE=/storage1/fs1/ha/Active/maherlab/hdang/mCRPC-omics/atac-seq/alignments
#Sample ID
ID=$1
#Working dir for sample
WD=/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/mCRPC/"$ID"

mkdir "$WD"
cd "$WD"


#FILTERING ALIGNMENTS 

#Remove duplicates from coordinate soorted BAMs
java -jar /usr/local/bin/picard.jar MarkDuplicates QUIET=true INPUT="$STORAGE"/"$ID".bowtie.sorted.bam OUTPUT="$WD"/"$ID".dedup.bam METRICS_FILE="$WD"/"$ID".dedup.metrics.txt VALIDATION_STRINGENCY=LENIENT \
TMP_DIR="$WD" REMOVE_DUPLICATES=true

#Remove secondary alignments, low MAPQ scores, mito reads, unmapped reads or mate unmapped reads, reads that fail platform, then sort by coord and index.
samtools view -h -q 30 -F 1804 -f 2 "$WD"/"$ID".dedup.bam | grep -v chrM | samtools sort -o "$WD"/"$ID".final.bam
#samtools index "$WD"/"$ID".final.bam

rm "$WD"/"$ID".dedup.bam

#Extracting insert/fragment sizes - column 9 of SAM for PE reads
#samtools view "$WD"/"$ID".final.bam | cut -f9 | sort | uniq -c | sort -b -k2,2n > "$WD"/"$ID".insertsize.txt

 



