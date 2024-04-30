##---------------------------------------------------------------------#
##
## Title: Train model on subset of SATB2+ SATB2- FACS sorted fetal data, and test on remaining saples
##
## Purpose of script: divide dataset into training and testing sets.
##                    1, Early/late train/ late test
##                    2. Early train / early test (x-fold validate)
##
##
##
## Author: Emma Walker
##
## Date Created: 26/04/2024
## updated from script 09/01/2024
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# 0. Define Parameters
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

# the top option only contains probes that are in both the adult and fetal datasets (~570000 probes)
# the bottom option is fetal only
#data <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/SFARI_MRC_merged_N_NN.rdat"
data <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/normalisedBetas_FACS_Fetal.rdat"

#nTrain <- 66

seed <- 1

early <- 10
mid <- 18



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

# source script to relax f stat p.value threshold from 1e-8 to 1e-7
source("scripts/Function_pickCompProbesMatrixRelaxed.r")



#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#


load(data)

SampleSheet$Basename <- as.character(SampleSheet$Basename)

# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


# subset to just satb2+ and satb-  fetal samples
SampleSheet <- SampleSheet[SampleSheet$Phenotype == "Fetal",]
SampleSheet <- SampleSheet[SampleSheet$Cell_Type %in% c("SATB2pos", "SATB2neg"),]

SampleSheet$Age <- as.numeric(as.character(SampleSheet$Age))

# add age bin col for plotting
SampleSheet$Age.bin <- rep(NA)
SampleSheet$Age.bin[SampleSheet$Age <= early] <- "Early"
SampleSheet$Age.bin[SampleSheet$Age > early & SampleSheet$Age <= mid] <- "Mid"
SampleSheet$Age.bin[SampleSheet$Age > mid] <- "Late"

table(SampleSheet$Age.bin)
# Early  Late   Mid - lifecourse data
#   19    15    45 

# Early  Late   Mid - fetal only data
#   19    19    45 

#----------------------------------------------------------------------#
# EARLY SAMPLES CROSS FOLD VALIDATION
#----------------------------------------------------------------------#

# Total of 19 samples
# Use 18 to train, and test 1
# loop over all samples

#set.seed(seed)

# split into early/mid/late
age <- "Late"

subPheno <- SampleSheet[SampleSheet$Age.bin == age,]
#subPheno <- subPheno[!subPheno$Basename == "203973690033_R04C01",]

subBetas <- as.matrix(betas[,subPheno$Basename])




modelOutput <- list()
for(i in 1:(nrow(subPheno)-1)){
  
  train <- subPheno[-c(i, i+1),]
  test <- subPheno[c(i, i+1),]
  
  trainBetas <- subBetas[,train$Basename] # check they are lined up
  # decrease to 1-5 threshold?
  # run rowFtests
  # plot p vals
  # subtract mean diff between satb2+ and -ve
  # are models selecting same sites in the models??
  
  # for early samples, calc mean across satb2+ satb2_ samplles and tabulate +ve against -ve
  #hyop, yper and hemi meth and same for later samples
  # for an indivdual - plot
  # how many switch hyper/hypo or other way round look at that against age - 
  # do per sample
  # are there more sites that could be cell type markers in the older vs younger
  # is hyper/hypo constant or erratic
  modelOutput[[i]] <- pickCompProbesMatrix(rawbetas = trainBetas,
                                     cellTypes = unique(train$Cell_Type),
                                     cellInd = train$Cell_Type,
                                     numProbes = 100,
                                     probeSelect = "auto")
  
  names(modelOutput)[[i]] <- paste0(test$Basename[1], "_", test$Basename[2])

}


# plot pvalues from model

allPvalPlots <- list()
counter <- 1
for(i in 1:length(modelOutput)){
  
  df <- modelOutput[[i]]
  compTab <- df$compTable
  
  compTab <- compTab[compTab$p.value < 0.00001,]
  
  p <- ggplot(compTab, aes(x=p.value))+
    geom_histogram()+
    ggtitle(paste0("nCoefsProbes =", sum(compTab$p.value < 1e-08)))+
    geom_vline(xintercept = 1e-08, colour="red")
    
   
  
  
  allPvalPlots[[counter]] <- p
  counter = counter + 1
}

pdf("plots/modelCeofPvalHistograms.pdf")
do.call(grid.arrange, c(allPvalPlots, ncol = 4))
dev.off()

######################################


# leave one out of each - then can create ulk tissues of differtn ratios 10,20,30

## not working for early samples... not getting any coef estimates!

## 

predPropAll <- list()
for(i in 1:length(modelOutput)){
  
  #test <- names(modelOutput)[[i]]
  testBetas <- as.matrix(subBetas[,test$Basename])
  model <- modelOutput[[i]]
  
  rInd<-rownames(testBetas)[rownames(testBetas) %in% rownames(model$coefEsts)]
  predPropAll[[i]] <- as.data.frame(projectCellTypeWithError(testBetas, model$coefEsts[rInd,]))
  
  names(predPropAll)[[i]] <- names(modelOutput)[[i]]
  
}


# not working as not getting a CETGYO score...


#----------------------------------------------------------------------#
# sanity check
#----------------------------------------------------------------------#
set.seed(1)
testTEST <- sample(SampleSheet$Basename, 2)
testTestpheno <- SampleSheet[SampleSheet$Basename %in% testTEST,]
testTestBetas <- as.matrix(betas[,testTEST])
row.names(testTestBetas) <- row.names(betas)

trainTest <- SampleSheet[!SampleSheet$Basename %in% testTEST,]
trainBetas <- as.matrix(betas[,colnames(betas) %in% trainTest$Basename])

testModel <- pickCompProbesMatrix(rawbetas = trainBetas,
                                  cellTypes = unique(trainTest$Cell_Type),
                                  cellInd = trainTest$Cell_Type,
                                  numProbes = 100,
                                  probeSelect = "any")

rInd<-rownames(testTestBetas)[rownames(testTestBetas) %in% rownames(testModel$coefEsts)]
testPredProp <- projectCellTypeWithError(testTestBetas, testModel$coefEsts[rInd,])



#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#


pCETYGO <- ggplot(predPropAll[[i]], aes(factor(0), CETYGO))+
  geom_boxplot()+
  ggtitle(project)+
  theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())








##########################################################################

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


