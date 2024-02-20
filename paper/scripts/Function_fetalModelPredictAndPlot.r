##---------------------------------------------------------------------#
##
## Title: Function for running fetal CETYGO models
##
## Purpose of script: Define function for running CETYGO fetal models on datasets
##
## Input: betas matrix, prjoect name (for titles), CETYGO model
## 
## Output: List object containing: model output matux, boxplot of CETYGO scores, boxplot of 
## CETYGO predictions
##
## Author: Emma Walker
##
## Date Created: 15/01/2024
##---------------------------------------------------------------------##


# define function to run predictor and output plots
fetalModelTest <- function(betas, project, model){
  
  ## identify which sites in the model are present in test data
  rInd<-rownames(betas)[rownames(betas) %in% rownames(model$coefEsts)]
  predPropCustom<-as.data.frame(projectCellTypeWithError(betas, model$coefEsts[rInd,]))
  
  
  # CETYGO score boxplot
  pCETYGO <- ggplot(predPropCustom, aes(factor(0), CETYGO))+
    geom_boxplot()+
    ggtitle(project)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
  
  
  # proportions boxplot
  plotdf <-as.data.frame(predPropCustom[,1:(ncol(predPropCustom)-2)])
  plotdf$Basename <- rownames(predPropCustom)
  plotdf <- reshape2::melt(plotdf)
  colnames(plotdf) <- c("Basename", "CellType", "Proportion")
  
  p <- ggplot(plotdf, aes(x=CellType, y = Proportion)) +
    geom_boxplot()+
    ggtitle(project)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  
  
  
  returnList <- list(predPropCustom, pCETYGO, p)
  
}
