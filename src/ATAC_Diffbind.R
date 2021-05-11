#!/bin/Rscript
#R script to analyze differential accessiblity for ATAC-seq data with Diffbind package
#docker(biowardrobe2/diffbind:v0.0.13)

#Load Diffbind library
library(DiffBind)

#Path for CSV doc to create diffbind object
doc <- "/storage1/fs1/ha/Active/maherlab/saha.d/projects/ATAC_seq/R_scripts/"

#Upload RDS file with samples and clinical annotations
samples <- readRDS(file = paste0(doc, "ATAC_mCRPC_DBA_samples.rds"))

#Create DBA object from csv sheet
print("Generating dba object from sample sheet")
atac_mcrpc <- dba(sampleSheet = samples)

saveRDS(object = atac_mcrpc, file = "ATAC_mCRPC_DBA_object.rds")

#Read in dba object with peaks/bams
atac_mcrpc <- readRDS(file = paste0(doc,"ATAC_mCRPC_DBA_object.rds" ))

#Generate count matrix with BAM files - can stop parallel processing if too much memory
print("Generating RPKM matrix")
atac_mcrpc <- dba.count(atac_mcrpc, bParallel = T, score = DBA_SCORE_RPKM, minOverlap = 0.5)

#Set contrasts and blocking variables (i.e covariates to regress out)
#print("Setting contrasts")
#atac_mcrpc <- dba.contrast(atac_mcrpc, categories = DBA_FACTOR, block = c(DBA_CONDITION, DBA_TISSUE))
#atac_mcrpc

#Perform DESeq2 testing with regularized GLM
#print("Starting Differential analysis")
#atac_mcrpc <-  dba.analyze(atac_mcrpc)

#Show results
#print("Differential analysis results")
#dba.show(atac_mcrpc, bContrasts = T)

#Save object to download
saveRDS(object = atac_mcrpc, file = "ATAC_mCRPC_DBA_results.rds")



