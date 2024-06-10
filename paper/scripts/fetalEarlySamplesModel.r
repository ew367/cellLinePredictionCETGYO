##---------------------------------------------------------------------#
##
## Title: early data leave one out
##
## Purpose of script: Create model using early data, exluding one satb2+/- sample
##                    at random
## 
##
## Author: Emma Walker
##
## Date Created: 05/06/2024
##---------------------------------------------------------------------##


#----------------------------------------------------------------------#
# 0. Define Parameters
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

data <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/normalisedBetas_FACS_Fetal.rdat"

early <- 10

flexThresh <- 1e-4

#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(CETYGO)
#library(factoextra)
#library(dplyr)
#library(stringr)

setwd(projDir)

# source script to relax f stat p.value threshold from 1e-8 to flexThresh
source("scripts/Function_pickCompProbesMatrixRelaxed.r")


#----------------------------------------------------------------------#
# IMPORT DATA
#----------------------------------------------------------------------#

load(data)

SampleSheet$Basename <- as.character(SampleSheet$Basename)
SampleSheet$Age <- as.numeric(as.character(SampleSheet$Age))


# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


# subset to just satb2+ and satb-  fetal samples
SampleSheet <- SampleSheet[SampleSheet$Phenotype == "Fetal",]
SampleSheet <- SampleSheet[SampleSheet$Cell_Type %in% c("SATB2pos", "SATB2neg"),]


# subset to just early samples
SampleSheet <- SampleSheet[SampleSheet$Age <= early,]
betas <- betas[,SampleSheet$Basename]


#----------------------------------------------------------------------#
# TRAIN MODEL
#----------------------------------------------------------------------#

# exclude one satb2- and one satb2= sample

# select test samples to exclude from model/to test on
set.seed(1)

test <- c(SampleSheet$Basename[sample(which(SampleSheet$Cell_Type == "SATB2pos"),1)],
          SampleSheet$Basename[sample(which(SampleSheet$Cell_Type == "SATB2neg"),1)])

testPheno <- SampleSheet[SampleSheet$Basename %in% test,]
testBetas <- betas[,testPheno$Basename]


trainPheno <- SampleSheet[!SampleSheet$Basename %in% test,]
trainBetas <- betas[,trainPheno$Basename]


#----------------------------------------------------------------------#
# CREATE MODEL
#----------------------------------------------------------------------#

modelOutput <- pickCompProbesMatrixRelaxed(rawbetas = trainBetas,
                                                cellTypes = unique(trainPheno$Cell_Type),
                                                cellInd = trainPheno$Cell_Type,
                                                numProbes = 100,
                                                probeSelect = "auto")


#----------------------------------------------------------------------#
# Run RowFtests on model betas - this is same output as modelOutput$compTable
#----------------------------------------------------------------------#

modelBetas <- trainBetas[row.names(modelOutput$coefEsts),]

rowFtests(modelBetas, as.factor(trainPheno$Cell_Type))


#----------------------------------------------------------------------#
# Mean diff between SATB2+ and SATB2-
#----------------------------------------------------------------------#

# all training sample probes

mean(rowMeans(betas[,trainPheno$Basename[trainPheno$Cell_Type == "SATB2pos"]]))
#0.5
mean(rowMeans(betas[,trainPheno$Basename[trainPheno$Cell_Type == "SATB2neg"]]))
#0.5


# just model probes
mean(rowMeans(modelBetas[,trainPheno$Basename[trainPheno$Cell_Type == "SATB2pos"]]))
#0.26
mean(rowMeans(modelBetas[,trainPheno$Basename[trainPheno$Cell_Type == "SATB2neg"]]))
#0.29


#----------------------------------------------------------------------#
# diff between SATB2+ and SATB2- from modelOutput compTable
#----------------------------------------------------------------------#

# all probes
diff <- abs(modelOutput$compTable$SATB2neg-modelOutput$compTable$SATB2pos)
range(diff)
# 3.327744e-08 3.358160e-01


# model probes
diff <- abs(modelOutput$compTable[row.names(modelOutput$coefEsts), "SATB2neg"] - 
            modelOutput$compTable[row.names(modelOutput$coefEsts), "SATB2pos"])
range(diff)
# [1] 0.007957241 0.197987740

# why is the probe with 0.3 difference not included in the model?



##notes

# for every satb2+ sample (satb2+ > satb2-)
# pair it in turn wih each satb2-ve sample
# this creates unique testing data each time
# and tests each sample 10(?) times
# can create bulk data of varying proportions w


