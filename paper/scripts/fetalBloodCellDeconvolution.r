##---------------------------------------------------------------------#
##
## Title: run blood deconvolution on bulk and fetal data
##
## Purpose of script: see whether the amount of blood cells changes with age
## in:
## 
## 1. fetal bulk (testing) data
## 2. fetal fans (training) data
##
## Author: Emma Walker
##
## Date Created: 04/04/2024
##
##---------------------------------------------------------------------##


#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(CETYGO)


#----------------------------------------------------------------------#
# BULK DATA
#----------------------------------------------------------------------#

# load bulk data
load("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/2_normalised/fetalBulk_EX3_23pcw_n91.rdat")

# run CETYGO with blood ref on bulk data
rowIndex<-rownames(betas)[rownames(betas) %in% rownames(modelBloodCoef)]
predProp<-as.data.frame(projectCellTypeWithError(betas, modelBloodCoef[rowIndex,]))


# plot by age

predProp$Basename <- row.names(predProp)

plotdf <- left_join(predProp, pheno %>% dplyr::select(Basename, PCW))
dfmelt <-reshape2::melt(plotdf %>% dplyr::select(-c(CETYGO, nCGmissing)), id.vars = c("Basename", "PCW"))

ggplot(plotdf, aes(x=PCW, y=CETYGO)) + 
  geom_point()+
  ggtitle("Bulk fetal - Blood deconvolution")

ggplot(dfmelt, aes(x=PCW, y=value, colour=variable)) + 
  geom_point()+
  ggtitle("Bulk fetal - Blood deconvolution")+
  geom_smooth(se=FALSE)+
  ylab("Proportion")


# remove everything from env

rm(list = ls())


#----------------------------------------------------------------------#
# FANS DATA
#----------------------------------------------------------------------#

# load FANS data
load("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/SFARI_MRC_merged_N_NN.rdat")

# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


# subset to just satb2+ and satb- fetal samples
fetal <- SampleSheet[SampleSheet$Phenotype == "Fetal",]
pheno <- fetal[fetal$Cell_Type %in% c("SATB2pos", "SATB2neg"),]

# extract betas for these samples
betas <- as.matrix(betas[,pheno$Basename])

# remove objects no longer needed
rm(list=setdiff(ls(), c("betas", "pheno")))



# run CETYGO with blood ref on bulk data
rowIndex<-rownames(betas)[rownames(betas) %in% rownames(modelBloodCoef)]
predProp<-as.data.frame(projectCellTypeWithError(betas, modelBloodCoef[rowIndex,]))


# plot by age

predProp$Basename <- row.names(predProp)

plotdf <- left_join(predProp, pheno %>% dplyr::select(Basename, Age))
dfmelt <-reshape2::melt(plotdf %>% dplyr::select(-c(CETYGO, nCGmissing)), id.vars = c("Basename", "Age"))

ggplot(plotdf, aes(x=Age, y=CETYGO)) + 
  geom_point()+
  ggtitle("FANS fetal - Blood deconvolution")

ggplot(dfmelt, aes(x=Age, y=value, colour=variable)) + 
  geom_point()+
  ggtitle("FANS fetal - Blood deconvolution")+
  geom_smooth(se=FALSE)+
  ylab("Proportion")





