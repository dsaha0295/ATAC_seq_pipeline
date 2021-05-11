#!/bin/bash
#Plot ATAC-seq sample enrichment at genomic regions using deeptools

WD=/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/Files/Final_bigwigs
HI=/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/Files/High_genes_TSShg38.bed
LO=/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/Files/Low_genes_TSShg38.bed
PCA=/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/Files/Prostate_GO_TSS_hg38.bed

#computeMatrix reference-point -S "$WD"/DTB-004-C.bw "$WD"/DTB-019-PRO-A.bw "$WD"/DTB-035-H.bw "$WD"/DTB-036-A.bw "$WD"/DTB-040-E.bw "$WD"/DTB-055-PRO-2C.bw "$WD"/DTB-059-A.bw "$WD"/DTB-060-E.bw "$WD"/DTB-063-B.bw "$WD"/DTB-071-A.bw "$WD"/DTB-090-PRO-B.bw "$WD"/DTB-111-PRO-A.bw "$WD"/DTB-127-A.bw "$WD"/DTB-128-B.bw "$WD"/DTB-129-B.bw "$WD"/DTB-138-D.bw "$WD"/DTB-167-A.bw "$WD"/DTB-167-PRO-B.bw "$WD"/DTB-172-A.bw "$WD"/DTB-176-A.bw "$WD"/DTB-183-A.bw "$WD"/DTB-222-F.bw "$WD"/DTB-265-PRO-B.bw "$WD"/PR-040-E.bw "$WD"/PR-081-A.bw \
#-R "$HI" "$LO" "$PCA" --referencePoint center -a 2000 -b 2000 -o mCRPC_ATAC_Hi_QC.gz 

#plotProfile -m mCRPC_ATAC_Hi_QC.gz -out mCRPC_ATAC_Hi_QC.png --plotType=fill --numPlotsPerRow 5

computeMatrix reference-point -S "$WD"/DTB-005-D.bw "$WD"/DTB-009-A.bw "$WD"/DTB-019-B.bw "$WD"/DTB-032-A.bw "$WD"/DTB-034-A.bw "$WD"/DTB-053-A.bw "$WD"/DTB-061-E.bw "$WD"/DTB-077-PRO-C.bw "$WD"/DTB-080-C.bw "$WD"/DTB-080-PRO-C.bw "$WD"/DTB-085-B.bw "$WD"/DTB-092-C.bw "$WD"/DTB-101-A.bw "$WD"/DTB-119-PRO-B.bw "$WD"/DTB-124-A.bw "$WD"/DTB-127-PRO-A.bw "$WD"/DTB-130-A.bw "$WD"/DTB-135-A.bw "$WD"/DTB-143-C.bw "$WD"/DTB-146-A.bw "$WD"/DTB-149-A.bw "$WD"/DTB-206-E.bw "$WD"/DTB-259-A.bw "$WD"/PR-056-A.bw "$WD"/PR-095-A.bw \
-R "$HI" "$LO" "$PCA" --referencePoint center -a 2000 -b 2000 -o mCRPC_ATAC_Med_QC.gz 

plotProfile -m mCRPC_ATAC_Med_QC.gz -out mCRPC_ATAC_Med_QC.png --plotType=fill --numPlotsPerRow 5

#computeMatrix reference-point -S "$WD"/DTB-008-D.bw "$WD"/DTB-011-B.bw "$WD"/DTB-022-A.bw "$WD"/DTB-030-C.bw "$WD"/DTB-037-A.bw "$WD"/DTB-044-D.bw "$WD"/DTB-064-E.bw "$WD"/DTB-069-E.bw "$WD"/DTB-074-D.bw "$WD"/DTB-083-A.bw "$WD"/DTB-091-B.bw "$WD"/DTB-097-C.bw "$WD"/DTB-097-PRO-A.bw "$WD"/DTB-110-C.bw "$WD"/DTB-125-A.bw "$WD"/DTB-126-A.bw "$WD"/DTB-140-A.bw "$WD"/DTB-141-A.bw "$WD"/DTB-148-A.bw "$WD"/DTB-156-B.bw "$WD"/DTB-170-A.bw "$WD"/DTB-173-A.bw "$WD"/DTB-176-PRO-A.bw "$WD"/DTB-218-A.bw "$WD"/PR-120-C.bw \
#-R "$HI" "$LO" "$PCA" --referencePoint center -a 2000 -b 2000 -o mCRPC_ATAC_Lo_QC.gz  


#plotProfile -m mCRPC_ATAC_Lo_QC.gz -out mCRPC_ATAC_Lo_QC.png --plotType=fill --numPlotsPerRow 5