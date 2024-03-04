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

bulkData <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/Bulk/2_normalised/fetalBulk_EX3_23pcw_n91.rdat"


#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

#library(dplyr)
library(ggplot2)
library(pheatmap)

setwd(projDir)

# load bulk fetal data
load(bulkData)

# add age bin col for plotting
pheno$Age.bin <- rep(NA)
pheno$Age.bin[pheno$PCW <= 10] <- "Early"
pheno$Age.bin[pheno$PCW >10 & pheno$PCW <=18] <- "Mid"
pheno$Age.bin[pheno$PCW > 18] <- "Late"

# change col types
pheno$M.median <- as.numeric(as.character(pheno$M.median))
pheno$U.median <- as.numeric(as.character(pheno$U.median))
pheno$Basename <- as.character(pheno$Basename)

# source function to predict and plot
source("scripts/Function_fetalModelPredictAndPlot.r")

# load models
load("models/allFetalSatb2Auto.rdat")
load("models/allFetalSatb2Any.rdat")
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

#### NOTE!!! this model no longer exists.. it has been replaced by:
# allFetalSatb2modelAuto

# run predictor trained on all available satb2 fetal samples
allSamplesModel <- fetalModelTest(betas, "AllSamplesModel", allFetalSatb2model)
pred <- allSamplesModel[[1]]

highCET <- pred[pred$CETYGO > 0.05,]

pheno$highCET <- rep(FALSE)
pheno$highCET[pheno$Basename %in% row.names(highCET)] <- TRUE 

# plot intensities again
ggplot(pheno, aes(x = M.median, y = U.median, colour = highCET))+
  geom_point()+
  xlab("Median M intensity")+
  ylab("Median U intensity")+
  ggtitle("Scatter plot of Signal Intensities coloured by high/low CETYGO")


#----------------------------------------------------------------------#
# Methylation density plots
#----------------------------------------------------------------------#

# extract beta values for probes in model
plotBetas <- betas[row.names(betas) %in% row.names(allFetalSatb2model[[1]]),]

densityPlot(plotBetas, main = "all fetal samples model probes", sampGroups = pheno$Age.bin)

# look at 67 probes that are in all 3 models
overlapProbes <- intersect(intersect(row.names(allFetalSatb2model[[1]]),
                                     row.names(weeks16to20Satb2Model[[1]])),
                                     row.names(weeks20Satb2Model[[1]]))

overlapBetas <- betas[row.names(betas) %in% overlapProbes,]

densityPlot(overlapBetas, main = "72 probes in all models", sampGroups = pheno$Age.bin)

#----------------------------------------------------------------------#
# heat maps
#----------------------------------------------------------------------#

colnames(plotBetas) <- paste0(colnames(plotBetas), "_", pheno$Age.bin)
#heatmap.2(plotBetas)

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

# add in whether they are in clust or not
pheno$outliers <- rep(FALSE)
pheno$outliers[pheno$Basename %in% tBasenames] <- TRUE

# plot densities
densityPlot(plotBetas, main = "probes in all models", sampGroups = pheno$outliers)


#----------------------------------------------------------------------#
# set up for heat maps using pheatmap 
# adapted from A. Franklins pheatmapExample.r
#----------------------------------------------------------------------#

plot_samp <- pheno
plot_samp$num <- 1:nrow(plot_samp)

# Create annotation matrix
nCuts <- nlevels(factor(plot_samp$Cell_Type))

annotation_col <- data.frame(CellType=plot_samp$Cell_Type, Age=plot_samp$PCW)
rownames(annotation_col) <- paste0(plot_samp$num, "_", plot_samp$Age.bin, "_", plot_samp$PCW)


#----------------------------------------------------------------------#
# pheatmap - auto probe select model
#----------------------------------------------------------------------#

# subset betas
plot_betas <- betas[row.names(betas) %in% row.names(allFetalSatb2modelAuto[[1]]),]
colnames(plot_betas) <- paste0(plot_samp$num, "_", plot_samp$Age.bin, "_", plot_samp$PCW)

# Run pheatmap
pheatmap(plot_betas, annotation_col=annotation_col, show_colnames=TRUE, show_rownames=FALSE, cutree_cols=nCuts, main="Bulk Fetal", angle_col=90,                           filename="plots/bulkAutoProbeSelectPheatmap.pdf", width=20, height=8)


#----------------------------------------------------------------------#
# pheatmap - any probe select model
#----------------------------------------------------------------------#

# subset betas
plot_betas <- betas[row.names(betas) %in% row.names(allFetalSatb2modelAny[[1]]),]
colnames(plot_betas) <- paste0(plot_samp$num, "_", plot_samp$Age.bin, "_", plot_samp$PCW)

# Run pheatmap
pheatmap(plot_betas, annotation_col=annotation_col, show_colnames=TRUE, show_rownames=FALSE, cutree_cols=nCuts, main="Bulk Fetal", angle_col=90,                           filename="plots/bulkAnyProbeSelectPheatmap.pdf", width=20, height=8)


#----------------------------------------------------------------------#
# box plots of CETYGO scores
#----------------------------------------------------------------------#
pred$Basename <- rownames(pred)
predBox <- left_join(pred, pheno %>% dplyr::select(Basename, outliers))

ggplot(predBox, aes(x=outliers, y=CETYGO)) + 
  geom_boxplot()


#----------------------------------------------------------------------#
# box plot of intensity ratios
#----------------------------------------------------------------------#
pheno$ratio <- as.numeric(as.character(pheno$ratio))

ggplot(pheno, aes(x=outliers, y=ratio)) + 
  geom_boxplot()

#no obvs differences 
