##---------------------------------------------------------------------#
##
## Title: DMRcate analysis
##
## Purpose of script: perform regional analysis of 
## schizophrenia vs controls EWAS results using DMRcate package
## Author: Emma Walker
##
## Date Created: 07-11-2023
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(DMRcate)
library(ENmix)


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
#args<-c("/lustre/projects/Research_Project-MRC190311/DNAm/MRC", "NeuN+")
dataDir <- args[1]
cellType <- args[2]

normData<-file.path(dataDir, "3_normalised/normalised.rdata")

outFile <- file.path(dataDir, paste0("CETYGO/EWAS/DMRcate/", cellType, "nullCTDMRs.rdat"))


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]
QCmetrics<-QCmetrics[!is.na(QCmetrics$CCDNAmAge),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# remove SNPs
celltypeNormbeta <- rmSNPandCH(celltypeNormbeta, dist=2, mafcut=0.05) #doesn't work


#----------------------------------------------------------------------#
# RUN DMRcate
#----------------------------------------------------------------------#

design <- model.matrix(~QCmetrics$Phenotype +  QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)

myAnnotation <- cpg.annotate(datatype = "array", celltypeNormbeta, analysis.type="differential", design=design, contrasts = FALSE, coef="QCmetrics$PhenotypeSchizophrenia")

DMRs = dmrcate(myAnnotation, lambda=1000, C=2)

results.ranges <- extractRanges(DMRs, genome = "hg19") 

save(results.ranges, file=outFile)
                                       
                                       
                                       
