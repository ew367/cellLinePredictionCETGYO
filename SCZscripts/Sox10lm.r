##---------------------------------------------------------------------#
##
## Title: EWAS with linear regression model for SOX10+ cell type
##
## Purpose of script: perform DNA methylation association analysis of 
## schizophrenia vs controls testing for main and cell-specific effects.
## seperately testing for cell type differences
##
## Author: Emma Walker (adapted from Eilis Hannon's script)
##
## Date Created: 2023-10-10
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# DEFINE ANALYSIS FUNCTION
#----------------------------------------------------------------------#

runEWAS<-function(row,QCmetrics){
  
  nullCT<-lm(row ~ QCmetrics$Phenotype +  QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
  
  # extract case control main effect
  return(c(summary(nullCT)$coefficients["QCmetrics$PhenotypeSchizophrenia", c(1,2,4)]))  
}



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
#args <- c("/lustre/projects/Research_Project-MRC190311/DNAm/MRC", "Sox10+")
dataDir <- args[1]
cellType <- args[2]


normData<-file.path(dataDir, "3_normalised/normalised.rdata")
CETYGOdata<-file.path(dataDir, "CETYGO/predictedProportions_MRC_SCZ.rdat")
resPath<-file.path(dataDir, "CETYGO/EWAS/")

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)
load(CETYGOdata)

print(paste0("running EWAS on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# take top 1000 rows for debugging
#betasSub <- celltypeNormbeta[1:10,]

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 3, byrow = TRUE)
#outtab<-matrix(data = parRapply(cl, betasSub, runEWAS, QCmetrics), ncol = 3, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
#rownames(outtab)<-rownames(betasSub)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P")  


save(outtab, file = file.path(paste0(resPath, cellType,"LM.rdata")))
