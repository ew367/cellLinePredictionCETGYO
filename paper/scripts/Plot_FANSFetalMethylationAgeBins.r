##---------------------------------------------------------------------#
##
## Title: Plot methylation profiles of fetal data coloured by age bins
##
## Purpose: plot bulk data intensities coloured by age bin
##        : plot bulk data intensities coloured by high/low CETYGO
##        : methylation density plots all samples model coloured by age bin
##        : methylation density plot of 72 probes in all models
##        : heatmaps of all sample model probes
##        : methylation density plot of samples that are clustering differently
## 
##
## Author: Emma Walker
##
## Date Created: 29/01/2024
##
##---------------------------------------------------------------------##

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

fansData <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/SFARI_MRC_merged_N_NN.rdat"


#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

#library(dplyr)
library(ggplot2)

setwd(projDir)

# load fans fetal data
load(fansData)

# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


# subset to just satb2+ and satb-  fetal samples
fetal <- SampleSheet[SampleSheet$Phenotype == "Fetal",]
satb2 <- fetal[fetal$Cell_Type %in% c("SATB2pos", "SATB2neg"),]
satb2$Basename <- as.character(satb2$Basename)


# load QC data

qc <- read.csv("samplesPassed_FACS_Fetal.csv", stringsAsFactors = F)
pheno <- left_join(satb2, qc %>% dplyr::select(Basename, M.median, U.median))

# match betas to pheno
betas <- betas[,pheno$Basename]


# add age bin col for plotting
pheno$Age.bin <- rep(NA)
pheno$Age.bin[pheno$Age <= 10] <- "Early"
pheno$Age.bin[pheno$Age >10 & pheno$Age <=18] <- "Mid"
pheno$Age.bin[pheno$Age > 18] <- "Late"


# source function to predict and plot
source("scripts/Function_fetalModelPredictAndPlot.r")

# load models
load("models/allFetalSatb2.rdat")
load("models/weeksTo20FetalSatb2.rdat")
load("models/weeks16to20FetalSatb2.rdat")


#----------------------------------------------------------------------#
# Plot intensity scores colourd by age bin
#----------------------------------------------------------------------#

# intensity plot
ggplot(pheno, aes(x = M.median, y = U.median, colour = Age.bin))+
  geom_point()+
  xlab("Median M intensity")+
  ylab("Median U intensity")+
  ggtitle("Scatter plot of Signal Intensities coloured by age bin")


#----------------------------------------------------------------------#
# Look at intensities of samples with high CETYGO score
#----------------------------------------------------------------------#


# run predictor trained on all available satb2 fetal samples
#allSamplesModel <- fetalModelTest(betas, "AllSamplesModel", allFetalSatb2model)
#pred <- allSamplesModel[[1]]

#highCET <- pred[pred$CETYGO > 0.05,]

#pheno$highCET <- rep(FALSE)
#pheno$highCET[pheno$Basename %in% row.names(highCET)] <- TRUE 

# plot intensities again
#ggplot(pheno, aes(x = M.median, y = U.median, colour = highCET))+
 # geom_point()+
  #xlab("Median M intensity")+
  #ylab("Median U intensity")+
  #ggtitle("Scatter plot of Signal Intensities coloured by high/low CETYGO")


#----------------------------------------------------------------------#
# Methylation density plots
#----------------------------------------------------------------------#

# extract beta values for probes in model
plotBetas <- as.matrix(betas[row.names(betas) %in% row.names(allFetalSatb2model[[1]]),])

densityPlot(plotBetas, main = "all fetal samples model probes", sampGroups = pheno$Age.bin)

# look at 67 probes that are in all 3 models

overlapProbes <- intersect(intersect(row.names(allFetalSatb2model[[1]]),
                                     row.names(weeks16to20Satb2Model[[1]])),
                                     row.names(weeks20Satb2Model[[1]]))

overlapBetas <- as.matrix(betas[row.names(betas) %in% overlapProbes,])

densityPlot(overlapBetas, main = "72 probes in all models", sampGroups = pheno$Age.bin)


#----------------------------------------------------------------------#
# heat maps
#----------------------------------------------------------------------#
plotBetas <- as.matrix(betas[row.names(betas) %in% row.names(allFetalSatb2model[[1]]),])

colnames(plotBetas) <- paste0(colnames(plotBetas), "_", pheno$Cell_Type)

heatmap.2(plotBetas)

# plot so can extract samples
pdf("plots/heatmap.pdf", width=30, height = 20)
row_clust <- hclust(dist(plotBetas, method = 'euclidean'), method = 'ward.D2')
out <- heatmap.2(
  plotBetas,
  Rowv = as.dendrogram(row_clust))
dev.off()

# get samples which are in different cluster
t <-colnames(plotBetas)[out$colInd][20:39]
tBasenames <- substring(t, 1, 19)


# check CETYGO scores of these samples

outliers <- pred[row.names(pred) %in% tBasenames,]

# box plot of CETYGO scores grouped by in cluster or not


# why are the methylation values of these so different? What is different about these samples?
# look at bscon/institute/plate/sex/ratio

# look at training data quality

# add in whether they are in clust or not
pheno$outliers <- rep(FALSE)
pheno$outliers[pheno$Basename %in% tBasenames] <- TRUE

# plot densities
densityPlot(plotBetas, main = "probes in all models", sampGroups = pheno$outliers)



#----------------------------------------------------------------------#
# box plots of CETYGO scores
#----------------------------------------------------------------------#
pred$Basename <- rownames(pred)
predBox <- left_join(pred, pheno %>% dplyr::select(Basename, outliers))

ggplot(predBox, aes(x=outliers, y=CETYGO)) + 
  geom_boxplot()

