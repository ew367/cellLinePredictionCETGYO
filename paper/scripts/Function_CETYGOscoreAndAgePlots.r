##---------------------------------------------------------------------#
##
## Title: Define function for plotting CETYGO output over age
##
## Author: Emma Walker
##
## Date Created: 12/01/2024
## updated from script 04/12/2023
##---------------------------------------------------------------------##



#----------------------------------------------------------------------#
# Define Function for CETYGO score/proportion ~ Age plots
#----------------------------------------------------------------------#

plotCETYGOage <- function(model, title){
  
  modelOutput <- model[[1]]
  modelOutput$Basename <- row.names(modelOutput)
  
  plotdf <- left_join(modelOutput, pheno %>% dplyr::select(Basename, PCW))
  dfmelt <-reshape2::melt(plotdf %>% dplyr::select(-c(CETYGO, nCGmissing)), id.vars = c("Basename", "PCW"))
  
  p1 <- ggplot(plotdf, aes(x=PCW, y=CETYGO)) + 
    geom_point()+
    geom_hline(yintercept=0.07, linetype="dashed", color = "red")+
    ggtitle(title)
  
  p2 <- ggplot(dfmelt, aes(x=PCW, y=value, colour=variable)) + 
    geom_point()+
    ggtitle(title)+
    geom_smooth(se=FALSE)+
    ylab("Proportion")+
    ylim(-2,2.1)+
    geom_hline(yintercept=0, linetype="dashed")+
    geom_hline(yintercept=1, linetype="dashed")
  
  returnList <- list(p1, p2)
}
