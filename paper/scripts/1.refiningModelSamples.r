##---------------------------------------------------------------------#
##
## Title: Model Samples Refinement
##
## Purpose of script: This script looks at various combinations of samples
##                    that could be included in the prediction models
##                    
##                    It uses PCA analysis and heat maps to visually inspect
##                    the data and decide on 3 three groups of samples to take
##                    forward for cross fold validation
##
##
##                The groups are: all samples
##                                13-28 week samples (large age range while
##                                having distinct satb2+/-)
##                                14-20 weeks (most seperation for satb2+ satb2-)  
##        
##                   
##
##---------------------------------------------------------------------#



#----------------------------------------------------------------------#
# Define Parameters
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

data <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/normalisedBetas_FACS_Fetal.rdat"

minAge <- 16
maxAge <- 20

toExclude <- ""


#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(CETYGO)
#library(factoextra)
#library(dplyr)
#library(stringr)

setwd(projDir)

manifest<-fread("/lustre/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv", skip=7, fill=TRUE, data.table=F)

#----------------------------------------------------------------------#
# IMPORT DATA
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

# remove objects no longer needed
#rm(list=setdiff(ls(), c("satb2", "satb2Betas")))

#remove sex specific probes
auto.probes<-manifest$IlmnID[manifest$CHR != "X" & manifest$CHR != "Y" & manifest$CHR != "MT"]

satb2Betas<-satb2Betas[row.names(satb2Betas) %in% auto.probes,]


#----------------------------------------------------------------------#
# Select Samples
#----------------------------------------------------------------------#


subPheno <- satb2[satb2$Age >= minAge & satb2$Age <= maxAge,]
subPheno <- subPheno[!subPheno$Individual_ID == 11947,]
subPheno <- subPheno[!subPheno$Basename == "203968030028_R06C01",]

subBetas <- as.matrix(satb2Betas[,subPheno$Basename])



#----------------------------------------------------------------------#
# PCA plots
#----------------------------------------------------------------------#


pca.res <- prcomp(t(subBetas))

plotdf <- as.data.frame(pca.res$x)
plotdf$Basename <- rownames(plotdf)
plotdf <- left_join(plotdf, subPheno)

var_explained <- pca.res$sdev^2/sum(pca.res$sdev^2)


# cell type
ggplot(plotdf, aes(x=PC1,y=PC2, colour = Cell_Type)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")

#Age
ggplot(plotdf, aes(x=PC1,y=PC2, colour = Age)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")

#Sex
ggplot(plotdf, aes(x=PC1,y=PC2, colour = Sex)) + geom_point(size=4) +
  labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
       y=paste0("PC2: ",round(var_explained[2]*100,1),"%")) +
  theme(legend.position="top")


#----------------------------------------------------------------------#
# PCA correlations
#----------------------------------------------------------------------#

pcaCorPlotDF <- cbind(plotdf[,1:5], weeks16to20 %>% dplyr::select(Sex, Age, Institute, Individual_ID, Plate_Location, Cell_Type, Tissue_Type))

pcaCorPlotDF$Sex <- as.numeric(pcaCorPlotDF$Sex)
pcaCorPlotDF$Institute <- as.numeric(as.factor(pcaCorPlotDF$Institute))
pcaCorPlotDF$Plate_Location <- as.numeric(as.factor(pcaCorPlotDF$Plate_Location))
pcaCorPlotDF$Individual_ID <- as.numeric(as.factor(pcaCorPlotDF$Individual_ID))
pcaCorPlotDF$Cell_Type <- as.numeric(as.factor(pcaCorPlotDF$Cell_Type))
pcaCorPlotDF$Tissue_Type <- as.numeric(as.factor(pcaCorPlotDF$Tissue_Type))

corrplot(cor(pcaCorPlotDF, use = "p"))


#----------------------------------------------------------------------#
# HEATMAPS
#----------------------------------------------------------------------#

# heat map
cellTypes<-sort(unique(subPheno$Cell_Type))

sample_anno<-subPheno[,c("Age","Sex", "Cell_Type")[c("Age","Sex", "Cell_Type") %in% colnames(subPheno)]]
if("Age" %in% colnames(sample_anno)){
  sample_anno$Age<-as.numeric(sample_anno$Age)
}
rownames(sample_anno)<-subPheno$Basename

sample_anno<-sample_anno[,colSums(!is.na(sample_anno))> 0, drop=F]

sigma<-apply(subBetas, 1, sd)



pheatmap(subBetas[order(sigma, decreasing = TRUE)[1:500],], annotation_col = sample_anno,  show_colnames = FALSE, show_rownames = FALSE, cutree_cols = length(cellTypes), main = "")

