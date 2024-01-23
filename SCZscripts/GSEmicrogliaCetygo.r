##---------------------------------------------------------------------#
##
## Title: Calculate CETYGO scores and estimatd cell proportions from normalised DNAm data
##
## Purpose of script: estimations of cell type proportions to be used as covariates in downstream analyses 
##
## Author: Emma Walker
##
## Date Created: November 2023
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# NOTES
#----------------------------------------------------------------------#
# project folder is provided on command line


#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#
args<-commandArgs(trailingOnly = TRUE)
#dataDir <- args[1]
dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO"

normData<-file.path(dataDir, "/data/GSE191200_files/QCSamples_Normalised.rdat")
resPath<-file.path(dataDir, "microglia")

#cells <- c("Double-", "NeuN+", "Sox10+")
project <- "microglia"

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("genefilter", "minfi"))
install.packages("quadprog")

#install devtools to install from GitHub
#install.packages("devtools")
library(devtools)
#install_github("ds420/CETYGO", ref = "add-brain-panels")
library(CETYGO)
library(gridExtra)


#----------------------------------------------------------------------#
# LOAD DATA
#----------------------------------------------------------------------#
setwd(dataDir)
load(normData)

# ensure sample sheet is in same order as data
identical(colnames(betas), as.character(SampleSheet$Basename))

# subset to just celltypes of interest
#QCmetricsCellTypes <- QCmetrics[which(QCmetrics$Cell.type %in% cells),]
#betasCellTypes <- celltypeNormbeta[,colnames(celltypeNormbeta) %in% QCmetricsCellTypes$Basename]

## see if any CETYGO data already exists
#if(file.exists(qcData)){
 # load(qcData)
  ## check contains all required samples
  #if(nrow(QCmetrics) == nrow(sampleSheet)){
   # print("QC file loaded")
  #} else {
   # QCmetrics<-sampleSheet
    #print("QC file to be updated with new samples")
  #}
#} else{
 # QCmetrics<-sampleSheet
  #print("QC object initiated")
#}




#----------------------------------------------------------------------#
# DEFINE FUNCTION TO RUN CETYGO
#----------------------------------------------------------------------#


adultBrainCETYGO <- function(betas, project){
  
  predPropAll<-list()
  CETYGOplots <- list()
  propPlots <- list()
  counter = 1
  maxCETYGO <- 0
  minCETYGO <- 1
  
  for(method in names(modelBrainCoef)){
    for(j in 1:length(modelBrainCoef[[method]])){
      if(!is.null(modelBrainCoef[[method]][[j]])){
        predPropAll[[method]][[j]]<-projectCellTypeWithError(betas, modelBrainCoef[[method]][[j]])
      }
    } 
  }
  
  
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[1]]
    for(j in 1:length(type)){
      plotdf <- as.data.frame(type[[j]])
      
      if (ncol(plotdf) > 0) {
        
        if (max(plotdf$CETYGO) > maxCETYGO){
          maxCETYGO <- max(plotdf$CETYGO)}
        
        if (min(plotdf$CETYGO) < minCETYGO){
          minCETYGO <- min(plotdf$CETYGO)}
      }
    }
  }
  
  for(i in 1:length(predPropAll)){
    type <- predPropAll[[1]]
    for(j in 1:length(type)){
      plotdf <- as.data.frame(type[[j]])
      
      if (ncol(plotdf) > 0) {
        
        # CETYGO score boxplot
        pCETYGO <- ggplot(plotdf, aes(factor(0), CETYGO))+
          geom_boxplot()+
          coord_cartesian(ylim=c(minCETYGO - 0.01, maxCETYGO + 0.01))+
          ggtitle(paste0(project, " - ", names(predPropAll[i]), "-", j))+
          theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(), plot.title = element_text(size=6))
        
        CETYGOplots[[counter]] <- pCETYGO
        
        # proportions boxplot
        plotdf <-as.data.frame(plotdf[,1:(ncol(plotdf)-2)])
        plotdf$Basename <- rownames(plotdf)
        plotdf <- reshape2::melt(plotdf)
        colnames(plotdf) <- c("Basename", "CellType", "Proportion")
        
        p <- ggplot(plotdf, aes(x=CellType, y = Proportion)) +
          geom_boxplot()+
          ggtitle(paste0(project, " - ", names(predPropAll[i]), "-", j))+
          theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
        
        propPlots[[counter]] <- p
        
        counter = counter + 1
      }
    }
  }
  
  CETYGOplotFile <- paste0(resPath, "/plots/CETYGO_", project, ".pdf")
  pdf(CETYGOplotFile)
  do.call(grid.arrange, c(CETYGOplots, ncol = 4))
  dev.off()
  
  propPlotsFile <- paste0(resPath, "/plots/cellProps_", project, ".pdf")
  pdf(propPlotsFile, height = 20, width = 20)
  do.call(grid.arrange, c(propPlots, ncol = 4))
  dev.off()
  
  predPropOutFile <- paste0(resPath, "/predictedProportions_", project, ".rdat")
  save(predPropAll, file = predPropOutFile)
  
}

adultBrainCETYGO(betas, "microglia")


# subset to each celltype of interest for plots
########## see lm script for extraction of data for plotting

load(paste0(resPath, "/predictedProportions_", project, ".rdat"))

#rm stuff no longer needed 

cells <- c("Double-", "NeuN+", "Sox10+")

## code as for loop??

QCdouble <- QCmetrics[which(QCmetrics$Cell.type == "Double-"),]
doubleBetas <- celltypeNormbeta[, colnames(celltypeNormbeta) %in% QCdouble$Basename]

adultBrainCETYGO(doubleBetas, "doubleNeg")


QCNeuN <- QCmetrics[which(QCmetrics$Cell.type == "NeuN+"),]
betasNeuN <- celltypeNormbeta[, colnames(celltypeNormbeta) %in% QCNeuN$Basename]

adultBrainCETYGO(betasNeuN, "NeuNpos")






