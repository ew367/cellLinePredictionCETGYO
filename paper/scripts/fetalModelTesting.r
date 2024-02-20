##---------------------------------------------------------------------#
##
## Title: Test SATB2+ SATB2- model on fetal data
##
## Purpose of script: Test both of the models trained in the fetalModelTraining.r script on the bulk ## data to see which one performe best (lowest CETYGO score and 'sensible' predictions)
## 
## allFetalSatb2Model is all the Sat2b+ and Satb2- fetal samples available (79 samples)
## weeks20Satb2Model is restricted to just samples with age < 20 weeks (70 samples)
## 
##
## Author: Emma Walker
##
## Date Created: 12/01/2024
## updated from script 04/12/2023
##---------------------------------------------------------------------##


#----------------------------------------------------------------------#
# 0. Define Parameters
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

bulkData <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/fetalBulk_EX3_23pcw_n91.rdat"

#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(dplyr)
library(ggplot2)

setwd(projDir)

# load models
load("models/allFetalSatb2.rdat")
load("models/weeksTo20FetalSatb2.rdat")
load("models/weeks16to20FetalSatb2.rdat")


# source function to predict and plot
source("scripts/Function_fetalModelPredictAndPlot.r")


#extract bulk brain samples
# load("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/EPICBrainLifecourse.rdat")
#bulk <- epic.pheno[which(epic.pheno$Study_Name == "Foetal Bulk Brain" & epic.pheno$AgeCat =="Prenatal"),]
#bulkBetas <- epic.betas[,bulk$Basename]

# load bulk fetal data
load(bulkData)

#----------------------------------------------------------------------#
# FETAL BULK DATA
#----------------------------------------------------------------------#


## run predictor trained on all available satb2 fetal samples
allSamplesModel <- fetalModelTest(betas, "AllSamplesModel", allFetalSatb2model)

pdf("plots/allSamplesModel.pdf")
print(allSamplesModel[[2]])
print(allSamplesModel[[3]])
dev.off()


## run predictor trained on satb2 fetal samples < 20 weeks
weeks20Model <- fetalModelTest(betas, "weeks20Model", weeks20Satb2Model)

pdf("plots/Weeks20Model.pdf")
print(weeks20Model[[2]])
print(weeks20Model[[3]])
dev.off() 


## run predictor trained on satb2 fetal samples >=16 & <= 20 weeks wiht no 11947
load("models/weeks16to20FetalSatb2No11947.rdat")
weeks16to20no11947Model <- fetalModelTest(betas, "weeks16to20Model", weeks16to20Satb2Model)

pdf("plots/Weeks16to20Model.pdf")
print(weeks16to20no11947Model[[2]])
print(weeks16to20no11947Model[[3]])
dev.off()



## run predictor trained on satb2 fetal samples >=16 & <= 20 weeks with no 
weeks16to20no11947Model <- fetalModelTest(betas, "weeks16to20no11947Model", weeks16to20Satb2)

pdf("plots/Weeks16to20Model.pdf")
print(weeks16to20Model[[2]])
print(weeks16to20Model[[3]])
dev.off()


#----------------------------------------------------------------------#
# CETYGO SCORE ~ AGE
#----------------------------------------------------------------------#

# plot CETYGO score against age

model <- weeks16to20Model
title <- "weeks 16 to 20Model"

modelOutput <- model[[1]]
modelOutput$Basename <- row.names(modelOutput)

plotdf <- left_join(modelOutput, pheno %>% dplyr::select(Basename, PCW))

ggplot(plotdf, aes(x=PCW, y=CETYGO)) + 
  geom_point()+
  ggtitle(title)



#----------------------------------------------------------------------#
# PROPORTIONS ~ AGE
#----------------------------------------------------------------------#

model <- weeks16to20Model
title <- "weeks 16 to 20Model"

modelOutput <- model[[1]]
modelOutput$Basename <- row.names(modelOutput)

plotdf <- left_join(modelOutput, pheno %>% dplyr::select(Basename, PCW))

dfmelt <-reshape2::melt(plotdf %>% dplyr::select(-c(CETYGO, nCGmissing)), id.vars = c("Basename", "PCW"))

ggplot(dfmelt, aes(x=PCW, y=value, colour=variable)) + 
  geom_point()+
  ggtitle(title)+
  geom_smooth(se=FALSE)+
  ylab("Proportion")
  




################################### other datasets


# DRI IPSC
load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/agedNeurons/2_normalised/Normalised_DRICambridge.rdat")

dri <- fetalModelTest(betas, "DRI")

pdf("plots/dri.pdf")
print(dri[[2]])
print(dri[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Sam THC - ShSy5

load("data/THC_rerun_normalised.rdat")

sam <- fetalModelTest(betas2, "THC-ShSy5")

pdf("plots/thc.pdf")
print(sam[[2]])
print(sam[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Jenny - neurons

load("data/DasenData_badprobesrm_Neuron.rdata")

jen <- fetalModelTest(data.3_N, "ipsc neurons")

pdf("plots/jennyNeurons.pdf")
print(jen[[2]])
print(jen[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# hypo - astrocytes

load("data/hypo_beta_pf_ds_hybremoved_bisremoved.RData")
  
astro <- fetalModelTest(betas, "hypo astrocytes")

pdf("plots/hypoAstro.pdf")
print(astro[[2]])
print(astro[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Jenny microglia

load("data/DasenData_badprobesrm_Microglia.rdata")

micro <- fetalModelTest(data.3_MG, "microglia")

pdf("plots/micro.pdf")
print(micro[[2]])
print(micro[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


#### Alice bulk fetal

betas <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/Fetal/Bulk/2_normalised/fetalBulk_EX3_23pcw_n91.csv")
row.names(betas) <- betas$X
betas <- betas[,-c(1)]

bulk <- fetalModelTest(betas, "bulk")
pdf("plots/fetalBulk.pdf")
print(bulk[[2]])
print(bulk[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


#### Alice epic v2 sorted ####### need to look at this more closely!!!

load("/lustre/projects/Research_Project-MRC190311/DNAm/SFARI2023/3_normalised/SFARI2023_betas.rdat")

for(i in unique(SampleSheet$Cell_type)){
cellBetas <- betas[,colnames(betas) %in% SampleSheet$Basename[which(SampleSheet$Cell_type == i)]]

out <- fetalModelTest(cellBetas, paste0("epicV2 sorted ", i))
print(out[[2]])
print(out[[3]])
}

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Leo data

# tolosa 

load("data/data.tolosa.Rdata")

tolosaIPSCbetas <- betas.tolosa[,colnames(betas.tolosa) %in% pheno.tolosa$Sample[which(pheno.tolosa$Cell_State == "iPSC")]]

tolosaIPSC <- fetalModelTest(tolosaIPSCbetas, "tolosa iPSC")
print(tolosaIPSC[[2]])
print(tolosaIPSC[[3]])


tolosaNbetas <- betas.tolosa[,colnames(betas.tolosa) %in% pheno.tolosa$Sample[which(pheno.tolosa$Cell_State == "Neuron")]]

tolosaN <- fetalModelTest(tolosaNbetas, "tolosa Neurons")
print(tolosaN[[2]])
print(tolosaN[[3]])

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# kiselev 

kiselevipscbetas <- betas.kiselev[,colnames(betas.kiselev) %in% pheno.kiselev$full_name[which(pheno.kiselev$Cell_state == "iPSC")]]

kiselevipsc <- fetalModelTest(kiselevipscbetas, "kiselev iPSC")
print(kiselevipsc[[2]])
print(kiselevipsc[[3]])



kiselevNbetas <- betas.kiselev[,colnames(betas.kiselev) %in% pheno.kiselev$full_name[which(pheno.kiselev$Cell_state == "Neuron")]]

kiselevN <- fetalModelTest(kiselevNbetas, "kiselev Neurons")
print(kiselevN[[2]])
print(kiselevN[[3]])

kiselevNbetas <- betas.kiselev[,colnames(betas.kiselev) %in% pheno.kiselev$full_name[which(pheno.kiselev$Cell_state == "Neuron")]]

kiselevN <- fetalModelTest(kiselevNbetas, "kiselev Neurons")
print(kiselevN[[2]])
print(kiselevN[[3]])


rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# kim 

##### hES ===

load("data/data.kim.Rdata")
pheno.kim$Cell_State <- as.character(pheno.kim$Cell_State)

kimHESbetas <- betas.kim[,colnames(betas.kim) %in% pheno.kim$Sample[which(pheno.kim$Cell_State == "hES")]]

kimhES <- fetalModelTest(kimHESbetas, "kim hES")
print(kimhES[[2]])
print(kimhES[[3]])


#Neuron

kimNbetas <- betas.kim[,colnames(betas.kim) %in% pheno.kim$Sample[which(pheno.kim$Cell_State == "Neuron")]]

kimN <- fetalModelTest(kimNbetas, "kim Neurons")
print(kimN[[2]])
print(kimN[[3]])


# NPC 


kimNPCbetas <- betas.kim[,colnames(betas.kim) %in% pheno.kim$Sample[which(pheno.kim$Cell_State == "NPC")]]

kimNPC <- fetalModelTest(kimNbetas, "kim NPC")
print(kimNPC[[2]])
print(kimNPC[[3]])


### TN samples from training datset

load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/Fetal/FACS/2_normalised/normalisedBetas_FACS_Fetal.rdat")

SampleSheet$Cell_Type <- as.character(SampleSheet$Cell_Type)
SampleSheet$Cell_Type <- trimws(SampleSheet$Cell_Type)


TNbetas <- betas[,colnames(betas) %in% SampleSheet$Basename[which(SampleSheet$Cell_Type == "TN")]]

TN <- fetalModelTest(TNbetas, "TN")
print(TN[[2]])
print(TN[[3]])


## Total Nuclei

TNbetas <- betas[,colnames(betas) %in% SampleSheet$Basename[which(SampleSheet$Cell_Type == "Total Nuclei")]]

TN <- fetalModelTest(TNbetas, "Total Nuclei")
print(TN[[2]])
print(TN[[3]])















