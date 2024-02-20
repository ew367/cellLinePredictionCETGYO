##---------------------------------------------------------------------#
##
## Title: Train model on subset of SATB2+ SATB2- FACS sorted fetal data, and test on remaining saples
##
## Purpose of script: divide dataset into training and tesitng set.
## Create models from all fetal FACS satb2 samples, and from just samples with age WPC < 20 weeks.
##
## Additionally look for patterns in the datasets (PCA clustering histograms etc.)
##
## Author: Emma Walker
##
## Date Created: 09/01/2024
## updated from script 19/05/2023
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# 0. Define Parameters
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

data <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/SFARI_MRC_merged_N_NN.rdat"

nTrain <- 66

seed <- 1


#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(CETYGO)
library(factoextra)
library(dplyr)
library(stringr)

setwd(projDir)

# source function to predict and plot
source("scripts/Function_fetalModelPredictAndPlot.r")

load(data)

# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


# subset to just satb2+ and satb-  fetal samples
fetal <- SampleSheet[SampleSheet$Phenotype == "Fetal",]
satb2 <- fetal[fetal$Cell_Type %in% c("SATB2pos", "SATB2neg"),]


#----------------------------------------------------------------------#
# SPLIT DATA
#----------------------------------------------------------------------#

set.seed(seed)

trainPheno <- satb2[sample(nrow(satb2), size=nTrain), ]
trainBetas <- as.matrix(betas[,trainPheno$Basename])

testPheno <- satb2[!satb2$Basename %in% trainPheno$Basename,]
testBetas <- as.matrix(betas[, testPheno$Basename])


#----------------------------------------------------------------------#
# TRAIN MODEL
#----------------------------------------------------------------------#

# Select the sites to form the basis of the deconvolution.
trainModel <- pickCompProbesMatrix(rawbetas = trainBetas,
                                           cellTypes = unique(trainPheno$Cell_Type),
                                           cellInd = trainPheno$Cell_Type,
                                           numProbes = 100,
                                           probeSelect = "auto")

save(trainModel, file="models/trainTestSubsetModel.rdat")



#----------------------------------------------------------------------#
# TEST MODEL
#----------------------------------------------------------------------#


## run predictor trained on all available satb2 fetal samples
testModel <- fetalModelTest(testBetas, "train/test subset Model", trainModel)

pdf("plots/allSamplesModel.pdf")
print(testModel[[2]])
print(testModel[[3]])
dev.off()


## check if predicting correctly

modelOutput <- testModel[[1]]
modelOutput$Basename <- row.names(modelOutput)

plotdf <- left_join(modelOutput, testPheno %>% dplyr::select(Basename, Age, Cell_Type))

plotdf$Prediction <- NA
plotdf$Prediction[plotdf$SATB2pos > plotdf$SATB2neg] <- "SATB2pos"
plotdf$Prediction[plotdf$SATB2pos < plotdf$SATB2neg] <- "SATB2neg"
plotdf$correctPred[plotdf$Cell_Type == plotdf$Prediction] <- TRUE
plotdf$correctPred[plotdf$Cell_Type != plotdf$Prediction] <- FALSE

## these samples are not working well: possibly because only 8/10/11 weeks
# 1.398515690	-0.38107733	0.07138077	0	203973690037_R04C01	8	SATB2pos	SATB2neg
# 0.584450363	0.39017109	0.04900585	0	205061420018_R06C01	10	SATB2pos	SATB2neg
# 0.297490834	0.70601991	0.05567645	0	205061420052_R03C01	11	SATB2neg	SATB2pos


dfmelt <-reshape2::melt(plotdf %>% dplyr::select(-c(CETYGO, nCGmissing, Prediction, Cell_Type)), id.vars = c("Basename", "Age", "correctPred"))

ggplot(dfmelt, aes(x=Age, y=value, colour = correctPred, shape=variable)) + 
  geom_point()
  ggtitle("Train/Test model")+
  ylab("Proportion")


