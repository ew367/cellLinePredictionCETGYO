##---------------------------------------------------------------------#
##
## Title: DMRcate analysis
##
## Purpose of script: perform regional analysis of 
## schizophrenia vs controls EWAS results using DMRcate package
## Author: Emma Walker
##
## NOTE: needs to be run on command line
##
## Date Created: 07-11-2023
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(DMRcate)
#library(ENmix)


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
#args<-c("/lustre/projects/Research_Project-MRC190311/DNAm/MRC", "NeuN+")
dataDir <- args[1]
cellType <- args[2]

CETYGOdata<-file.path(dataDir, "CETYGO/predictedProportions_MRC_SCZ.rdat")
normData<-file.path(dataDir, "3_normalised/normalised.rdata")

outFile <- file.path(dataDir, paste0("CETYGO/EWAS/DMRcate/", cellType, "modelLMDMRs.rdat"))

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)
load(CETYGOdata)

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]
QCmetrics<-QCmetrics[!is.na(QCmetrics$CCDNAmAge),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# add CETYGO proportions
if (cellType == "Double-"){
  QCmetrics$Cell.Proportions <- predPropAll$IDOL[[5]][QCmetrics$Basename, "NeuNNeg_Sox10Neg_IRF8Pos"]
} else if (cellType == "NeuN+"){
  QCmetrics$Cell.Proportions <- predPropAll$ANOVA[[6]][QCmetrics$Basename, "NeuNPos_SOX6Pos"]
} else if (cellType == "Sox10+"){
  stop("Run separate script for Sox10+ data")
} else {
  stop("Cell Type not recognised....")
}


# remove SNPs
celltypeNormbeta <- rmSNPandCH(celltypeNormbeta, dist=2, mafcut=0.05)


#----------------------------------------------------------------------#
# RUN DMRcate
#----------------------------------------------------------------------#

design <- model.matrix(~QCmetrics$Phenotype + QCmetrics$Cell.Proportions + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)

myAnnotation <- cpg.annotate(datatype = "array", celltypeNormbeta, analysis.type="differential", design=design, contrasts = FALSE, coef="QCmetrics$PhenotypeSchizophrenia")

DMRs = dmrcate(myAnnotation, lambda=1000, C=2)

results.ranges <- extractRanges(DMRs, genome = "hg19") 

save(results.ranges, file=outFile)

#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#
groups <- c(Control="magenta", Schizophrenia="forestgreen")
cols <- groups[as.character(QCmetrics$Phenotype)]

#pdf("CETYGO/plots/DMRcateNeuNTestPlots.pdf", width=60, height=100)
#par(mfrow=c(1,1))
#DMR.plot(ranges=results.ranges, dmr=1, CpGs=celltypeNormbeta, phen.col=cols, what="Beta",
         #arraytype="EPIC", genome="hg19", toscale=TRUE, plotmedians=TRUE)
#dev.off()

