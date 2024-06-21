##---------------------------------------------------------------------#
##
## Title: FANS data PCA plots
##
## Purpose of script: Plot PC1 against Age for the FANS fetal data to
##                    visualise where satb2+ and satb2- differentiate
##
##                                          
##                    Use PC1 values to define cut offs for inclusion
##                    in model
##
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# Some of the code in this script (including pcaFunctions.r) has been 
# adapted from a script written by AF for the fetal lifecourse project


#----------------------------------------------------------------------#
# Define Parameters
#----------------------------------------------------------------------#

library(ggplot2)
library(data.table)

# load PCA scripts
source("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/pcaFunctions.r")


projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

data <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/normalisedBetas_FACS_Fetal.rdat"

posThresh <- 25  # needs to be < this value
negThresh <- 0 # needs to be < this value


#----------------------------------------------------------------------#
# Import Data
#----------------------------------------------------------------------#

load(data)

# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


SampleSheet$Basename <- as.character(SampleSheet$Basename)
SampleSheet$Age <- as.numeric(as.character(SampleSheet$Age))


# subset to just satb2+ and satb- fetal samples
fetal <- SampleSheet[SampleSheet$Phenotype == "Fetal",]
satb2 <- fetal[fetal$Cell_Type %in% c("SATB2pos", "SATB2neg"),]

# extract betas for these samples
satb2Betas <- as.matrix(betas[,satb2$Basename])


betas <- satb2Betas
pheno <- satb2
pheno$Age <- as.numeric(pheno$Age)


#remove sex probes
probes <- data.frame(Probe=rownames(betas))
rownames(probes) <- rownames(betas)
epicManifest<-fread("/lustre/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B4.csv", skip=7, fill=TRUE, data.table=F)
epicMan<-epicManifest[match(rownames(probes), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
probes <- cbind(probes, as.data.frame(epicMan))
XY <- which(probes$CHR %in% c('X','Y'))
betas.noXY <- betas[-XY,]# [1] 686929     83



#----------------------------------------------------------------------#
# Run PCA and plot
#----------------------------------------------------------------------#


# run PCA - all samples
pca.all <- pca.Gene(pheno, betas.noXY, varProbes=TRUE, nVarprobes=10000, scale=TRUE)

# combine pheno and PCA output into dataframe for plotting
plotdf <- cbind(pheno, pca.all$x)

# plot PC1 against Age, coloured by cell type
ggplot(plotdf, aes(x=Age, y=PC1, colour = Cell_Type))+
  geom_point()+
  geom_smooth()+
  scale_x_continuous(breaks=seq(0,44,2))+
  geom_vline(xintercept = 12.5, linetype = "dashed")
  


#----------------------------------------------------------------------#
# Define group of sample based on PC1 values for using in cross fold
#  validation and replot
#----------------------------------------------------------------------#

# remove samples that are beyond threshold vals
pass <- c(plotdf$Basename[which(plotdf$Cell_Type == "SATB2pos" & plotdf$PC1 > posThresh)],
          plotdf$Basename[which(plotdf$Cell_Type == "SATB2neg" & plotdf$PC1 < negThresh)])
          
#plotdfPass[,c("Cell_Type", "PC1", "Age")]


# subset to 'passed' samples only
passPheno <- pheno[pheno$Basename %in% pass,]
passBetas <- betas.noXY[,passPheno$Basename]

# re run PCA on passed samples only
pca.Pass <- pca.Gene(passPheno, passBetas, varProbes=TRUE, nVarprobes=10000, scale=TRUE)

# combine pheno and PCA output into dataframe for plotting
plotdfPass <- cbind(passPheno, pca.Pass$x)

# plot PC1 against Age, coloured by cell type
ggplot(plotdfPass, aes(x=Age, y=PC1, colour = Cell_Type))+
  geom_point()+
  geom_smooth()+
  scale_x_continuous(breaks=seq(0,44,2))+
  geom_vline(xintercept = 12.5, linetype = "dashed")



#----------------------------------------------------------------------#
# Save for downstream analysis
#----------------------------------------------------------------------#

save(passPheno, passBetas, file = "PCADefinedSamples.rdat")
