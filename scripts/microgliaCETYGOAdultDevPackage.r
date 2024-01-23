#Emma Walker
#22/11/2023
#E.M.Walker@exeter.ac.uk

## on command line "R version 4.2.1 (2022-06-23)"
# not working on Rstudio "R version 3.6.0 (2019-04-26)"
# or R version 3.6 on command line R

#### NOT THI SCRIPT!!! GSEmicrogliaCetygo is the one


if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("genefilter", "minfi"))
install.packages("quadprog")

#install devtools to install from GitHub
install.packages("devtools")
library(devtools)
install_github("ds420/CETYGO", ref = "add-brain-panels")
library(CETYGO)

setwd("/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/")

##### idol 5 is best model for microglia

# https://github.com/ds420/CETYGO/wiki/Estimating-cellular-composition-of-bulk-brain-tissue---development


# load test data

load("data/GSE191200_files/QC/Samples_Normalised.rdat")


# run on all models etc

predPropAll<-list()

for(method in names(modelBrainCoef)){
  for(j in 1:length(modelBrainCoef[[method]])){
    if(!is.null(modelBrainCoef[[method]][[j]])){
      predPropAll[[method]][[j]]<-projectCellTypeWithError(betas, modelBrainCoef[[method]][[j]])
    }
  } 
}

save(predPropAll, file = "CETYGOAdultBrainOutput.rdat")

pdf("plots/microgliaCETYGOboxplots.pdf")
for(i in 1:length(predPropAll)){
  type <- predPropAll[[1]]
  for(j in 1:length(type)){
    plotdf <- as.data.frame(type[[j]])
    if (ncol(plotdf) > 0) {
  boxplot(plotdf$CETYGO, main = paste0(names(predPropAll[i]), "-", j))
  }
  }
}

dev.off()
