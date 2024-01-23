##---------------------------------------------------------------------#
##
## Title: EWAS with linear regression model for individual cell types
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
  
  modelLM<-lm(row ~ QCmetrics$Phenotype + QCmetrics$Cell.Proportions + QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
  nullCT<-lm(row ~ QCmetrics$Phenotype +  QCmetrics$CCDNAmAge + QCmetrics$Sex + QCmetrics$Tissue.Centre)
  
  # extract case control main effect and cell proportion effect
  return(c(summary(modelLM)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)],
           summary(modelLM)$coefficients["QCmetrics$Cell.Proportions",c(1,2,4)],
           
           # extract cell specific case control effect
           summary(nullCT)$coefficients["QCmetrics$PhenotypeSchizophrenia",c(1,2,4)]))
}



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
#dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/MRC"
cellType <- args[2]
#cellType <- "NeuN+"

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
    
 

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]

# take top 100 rows for debugging
#betasSub <- celltypeNormbeta[1:100,]

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)
clusterExport(cl, list("runEWAS"))

outtab<-matrix(data = parRapply(cl, celltypeNormbeta, runEWAS, QCmetrics), ncol = 9, byrow = TRUE)
#outtab<-matrix(data = parRapply(cl, betasSub, runEWAS, QCmetrics), ncol = 9, byrow = TRUE)


rownames(outtab)<-rownames(celltypeNormbeta)
#rownames(outtab)<-rownames(betasSub)
colnames(outtab)<-c("SCZ_coeff", "SCZ_SE", "SCZ_P", 
                    paste0(cellType,"_coeff"), paste0(cellType,"_SE"), paste0(cellType,"_P"),
                    "nullCT_SCZ_coeff", "nullCT_SCZ_SE", "nullCT_SCZ_P") 


save(outtab, file = file.path(paste0(resPath, cellType,"LM.rdata")))
