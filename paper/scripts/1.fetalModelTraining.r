##---------------------------------------------------------------------#
##
## Title: Train model on SATB2+ SATB2- FACS sorted fetal data
##
## Purpose of script: Create models from:
##
## 1. all fetal satb2 samples using probe select "any" and "both"
## 2. just samples with age WPC < 20 weeks
## 3. samples age > 16 weeks and < 20 weeks
## 4. excluding outlier sample from above model
##
## Additionally look for patterns in the datasets:
## - PCA clustering
## - histograms
##
## Author: Emma Walker
##
## Date Created: 09/01/2024
## updated from script 19/05/2023
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# github CETYGO tutorial:

# https://github.com/ds420/CETYGO/blob/main/vignettes/QuantifyErrorInCellularCompositionEstimate.Rmd


# minfi::estimateCellCounts

# https://rdrr.io/bioc/minfi/man/estimateCellCounts.html

#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(CETYGO)
library(factoextra)
library(dplyr)
library(stringr)

setwd("/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/")

# load in fetal FANs data
# this is the file path from gpfs used in the previous version of the script, which is currently offline
#load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/Fetal/FACS/2_normalised/normalisedBetas_FACS_Fetal.rdat")


#load("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/SFARI_MRC_merged_N_NN.rdat")

load("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/SFARI_MRC_merged_N_NN.rdat")


# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


# subset to just satb2+ and satb- fetal samples
fetal <- SampleSheet[SampleSheet$Phenotype == "Fetal",]
satb2 <- fetal[fetal$Cell_Type %in% c("SATB2pos", "SATB2neg"),]

# extract betas for these samples
satb2Betas <- as.matrix(betas[,satb2$Basename])

# remove objects no longer needed
rm(list=setdiff(ls(), c("satb2", "satb2Betas")))



#----------------------------------------------------------------------#
# TRAIN MODEL ON ALL SATB2 SAMPLES 
#----------------------------------------------------------------------#

# using probeSelect="auto"
allFetalSatb2modelAuto <- pickCompProbesMatrix(rawbetas = satb2Betas,
                                           cellTypes = unique(satb2$Cell_Type),
                                           cellInd = satb2$Cell_Type,
                                           numProbes = 100,
                                           probeSelect = "auto")

save(allFetalSatb2modelAuto, file="models/allFetalSatb2Auto.rdat")



#using probeSelect="any"
allFetalSatb2modelAny <- pickCompProbesMatrix(rawbetas = satb2Betas,
                                           cellTypes = unique(satb2$Cell_Type),
                                           cellInd = satb2$Cell_Type,
                                           numProbes = 100,
                                           probeSelect = "any")

save(allFetalSatb2modelAny, file="models/allFetalSatb2Any.rdat")


#compare probe lists

auto <- rownames(allFetalSatb2modelAuto[[1]])
any <- rownames(allFetalSatb2modelAny[[1]])

inBoth <- intersect(auto, any)
# 100 - presumably hypo/hyper?


# plot densities

densityPlot(satb2Betas[rownames(satb2Betas) %in% auto,])

densityPlot(satb2Betas[rownames(satb2Betas) %in% any,])



#----------------------------------------------------------------------#
# TRAIN MODEL ON SAMPLES SAMPLES < 20 weeks
#----------------------------------------------------------------------#

# subset to samples < 20 weeks
weeks20 <- satb2[satb2$Age < 20,]
weeks20betas <- as.matrix(betas[,weeks20$Basename])


# Select the sites to form the basis of the deconvolution.

weeks20Satb2Model <- pickCompProbesMatrix(rawbetas = weeks20betas,
                                          cellTypes = unique(weeks20$Cell_Type),
                                          cellInd = weeks20$Cell_Type,
                                          numProbes = 100,
                                          probeSelect = "auto")

save(weeks20Satb2Model, file="models/weeksTo20FetalSatb2.rdat")



#----------------------------------------------------------------------#
# TRAIN MODEL ON SAMPLES SAMPLES > 16 weeks and < 20 weeks
#----------------------------------------------------------------------#

## may also want to exclude 11947?


# subset to samples < 20 weeks
weeks16to20 <- satb2[satb2$Age >= 16 & satb2$Age <= 20,]
weeks16to20betas <- as.matrix(betas[,weeks16to20$Basename])


# Select the sites to form the basis of the deconvolution.

weeks16to20Satb2Model <- pickCompProbesMatrix(rawbetas = weeks16to20betas,
                                          cellTypes = unique(weeks16to20$Cell_Type),
                                          cellInd = weeks16to20$Cell_Type,
                                          numProbes = 100,
                                          probeSelect = "auto")

save(weeks16to20Satb2Model, file="models/weeks16to20FetalSatb2.rdat")



#----------------------------------------------------------------------#
# EXCLUDE 19947 FROM ABOVE MODEL
#----------------------------------------------------------------------#


weeks16to20 <- weeks16to20[!weeks16to20$Individual_ID == 11947,]
weeks16to20betas <- as.matrix(betas[,weeks16to20$Basename])


# Select the sites to form the basis of the deconvolution.
weeks16to20Satb2Model <- pickCompProbesMatrix(rawbetas = weeks16to20betas,
                                              cellTypes = unique(weeks16to20$Cell_Type),
                                              cellInd = weeks16to20$Cell_Type,
                                              numProbes = 100,
                                              probeSelect = "auto")

save(weeks16to20Satb2Model, file="models/weeks16to20FetalSatb2No11947.rdat")




#----------------------------------------------------------------------#
# SOME DATA EXPLORATION
#----------------------------------------------------------------------#


hist(weeks20$Age)

# make sure sample order matches in pheno and betas
print(identical(as.character(weeks20$Basename), colnames(weeks20betas)))
#[1] TRUE



### look how the data clusters

pca.res <- prcomp(t(weeks20betas))

#scree plot
fviz_eig(pca.res)

var_explained <- pca.res$sdev^2/sum(pca.res$sdev^2)

plotdf <- as.data.frame(pca.res$x)
plotdf$Basename <- rownames(plotdf)
plotdf <- left_join(plotdf, weeks20)



# cell type
ggplot(plotdf, aes(x=PC1,y=PC2, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")


ggplot(plotdf, aes(x=PC1,y=PC3, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(legend.position="top")



ggplot(plotdf, aes(x=PC1,y=PC4, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC4: ",round(var_explained[4]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC1,y=PC5, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC5: ",round(var_explained[5]*100,1),"%")) +
  theme(legend.position="top")


ggplot(plotdf, aes(x=PC2,y=PC3, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(legend.position="top")


ggplot(plotdf, aes(x=PC2,y=PC4, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC4: ",round(var_explained[4]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC2,y=PC5, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC5: ",round(var_explained[5]*100,1),"%")) +
  theme(legend.position="top")


## PC 3 and 4 seperate sex

ggplot(plotdf, aes(x=PC3,y=PC4, colour = Sex)) + geom_point(size=4) +
  labs(x=paste0("PC3: ",round(var_explained[3]*100,1),"%"),
       y=paste0("PC4: ",round(var_explained[4]*100,1),"%")) +
  theme(legend.position="top")


ggplot(plotdf, aes(x=PC3,y=PC5, colour = Sex)) + geom_point(size=4) +
  labs(x=paste0("PC3: ",round(var_explained[3]*100,1),"%"),
       y=paste0("PC5: ",round(var_explained[5]*100,1),"%")) +
  theme(legend.position="top")



# hist of age

hist(satb2$Age)

#----------------------------------------------------------------------#
# SAMPLES  > 16 weeks & < 20 weeks
#----------------------------------------------------------------------#

# check cluctering on mid range sample
nrow(satb2[satb2$Age >= 16 & satb2$Age <= 20,])

mid <- satb2[satb2$Age >= 16 & satb2$Age <= 20,]
# exclude sample 11947 
mid <- mid[!mid$Individual_ID == 11947,]
midbetas <- betas[,mid$Basename]

# make sure sample order matches in pheno and betas
print(identical(as.character(mid$Basename), colnames(midbetas)))
#[1] TRUE



### look how the data clusters

pca.res <- prcomp(t(midbetas))

#scree plot
fviz_eig(pca.res)

var_explained <- pca.res$sdev^2/sum(pca.res$sdev^2)

plotdf <- as.data.frame(pca.res$x)
plotdf$Basename <- rownames(plotdf)
plotdf <- left_join(plotdf, mid)



# cell type
ggplot(plotdf, aes(x=PC1,y=PC2, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC1,y=PC3, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(legend.position="top")



ggplot(plotdf, aes(x=PC2,y=PC3, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC2,y=PC4, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC4: ",round(var_explained[4]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC2,y=PC5, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC5: ",round(var_explained[5]*100,1),"%")) +
  theme(legend.position="top")


#----------------------------------------------------------------------#
# SAMPLES < 16 weeks
#----------------------------------------------------------------------#

# check cluctering on mid range sample
nrow(satb2[satb2$Age <= 16,])

weeks16 <- satb2[satb2$Age <= 16,]
# exclude sample 11947 
#mid <- mid[!mid$Individual_ID == 11947,]
weeks16betas <- betas[,weeks16$Basename]

pca.res <- prcomp(t(weeks16betas))

#scree plot
fviz_eig(pca.res)

var_explained <- pca.res$sdev^2/sum(pca.res$sdev^2)

plotdf <- as.data.frame(pca.res$x)
plotdf$Basename <- rownames(plotdf)
plotdf <- left_join(plotdf, weeks16)



# cell type
ggplot(plotdf, aes(x=PC1,y=PC2, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC1,y=PC3, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC1,y=PC4, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC4: ",round(var_explained[4]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC2,y=PC3, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC3: ",round(var_explained[3]*100,1),"%")) +
  theme(legend.position="top")

ggplot(plotdf, aes(x=PC2,y=PC4, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC2: ",round(var_explained[2]*100,1),"%"),
       y=paste0("PC4: ",round(var_explained[4]*100,1),"%")) +
  theme(legend.position="top")
