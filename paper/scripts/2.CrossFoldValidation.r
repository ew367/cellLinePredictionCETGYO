##---------------------------------------------------------------------#
##
## Title: Train models
##
## Purpose of script: Train CETYGO models using the 3 groups of samples
##                    chosen from script 1.refining model samples and..
##                    Perform cross fold validation
##                    
##                    The groups are: all samples
##                                    13-28 week samples (large age range while
##                                      having distinct satb2+/-)
##                                    14-20 weeks (most seperation for satb2+
##                                      satb2-)  
##                                    PCA samples from script 1.b
##        
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# only include matched samples? - check how many are matched/unmatched
# need a way to specifify PCA samples - 
# 3rd arg - PCA = TRUE/FALSE - if true load data from 1.b

# autosomes only

#----------------------------------------------------------------------#
# Define Parameters
#----------------------------------------------------------------------#

#read in args here and use to run on the different samples

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"

data <- "/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/2_normalised/normalisedBetas_FACS_Fetal.rdat"

manifest<-fread("/lustre/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv", skip=7, fill=TRUE, data.table=F)

#minAge <- args[1]
minAge <- 14

#maxAge <- args[2]
maxAge <- 20

# this only wants to be for not all samples models
if(minAge != 0){
  toExclude <- c("203968030028_R06C01", "11947")
}


# if PCA=TRUE

#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(CETYGO)

setwd(projDir)


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


# remove objects no longer needed
#rm(list=setdiff(ls(), c("satb2", "satb2Betas")))

#----------------------------------------------------------------------#
# Select samples to leave out in cross fold validation
#----------------------------------------------------------------------##

# n.b. some samples appear many times

# take individuals that have one satb2+ and one satb2- sample
leaveOut <- c() 
for(i in unique(subPheno$Individual_ID)){
  indivCT <- sort(subPheno$Cell_Type[subPheno$Individual_ID == i], decreasing = F)
  if(length(indivCT) == 2){
    if(identical(c("SATB2neg","SATB2pos"), indivCT)){
      leaveOut <- c(leaveOut, i)
    }
  }
}


# randomly select 10 individuals to leave out
set.seed(1)
leaveOut <- sample(leaveOut,10)


#----------------------------------------------------------------------#
# Create model and run cross fold validation
#----------------------------------------------------------------------#

# for each in leaveOut
# create train and test data
# create model using train data
# use test data to create pseudo bulk data
# run model on pseudo bulk data

i="13359"

modelOutput <- list()
for(i in leaveOut){
  
  trainPheno <- subPheno[!subPheno$Individual_ID == i,]
  trainBetas <- betas[,trainPheno$Basename]
  
  modelOutput[[i]] <- pickCompProbesMatrix(rawbetas = trainBetas,
                                          cellTypes = unique(trainPheno$Cell_Type),
                                            cellInd = trainPheno$Cell_Type,
                                          numProbes = 100,
                                        probeSelect = "any")
  
}


# for 1 model/individual left out combo
# note - don't need to do 100% pos as already have this
# rename 100% neg sample for

#extract probes used in training model
modelProbes <- row.names(modelOutput[[1]]$coefEsts)

#get testing data pheno info
testPheno <- subPheno[which(subPheno$Individual_ID == i),]

#get testing data betas and order so satbpos is 1st col, and neg is 2nd
testProbes <- betas[modelProbes,
                    c(testPheno$Basename[testPheno$Cell_Type == "SATB2pos"],
                      testPheno$Basename[testPheno$Cell_Type == "SATB2neg"])]

 
#create pseudoBulk samples - (10% pos + 90% neg, 20% pos + 80% neg, ...)
pseudoBulk <- testProbes[,2] # satb2-ve sample

for(j in 1:10){
  
  pBulk <- testProbes[,1]/10*j + testProbes[,2]/10*(10-j)
  pseudoBulk <- cbind(pseudoBulk, pBulk)
  
}

# colnames are % satb2+
for(k in 1:10){
  #print(paste0(colnames(testProbes)[k+1], k*10))
  colnames(pseudoBulk)[k+1] <- paste0(colnames(pseudoBulk)[k+1], k*10)
}
colnames(pseudoBulk)[1] <- "pBulk0"

# run predictor on pseudo bulk samples

#rInd<-rownames(betas)[rownames(betas) %in% rownames(model$coefEsts)]
#predPropCustom<-as.data.frame(projectCellTypeWithError(betas, model$coefEsts[rInd,]))

predPropCustom<-as.data.frame(projectCellTypeWithError(pseudoBulk, modelOutput[[1]]$coefEsts))


#----------------------------------------------------------------------#
# Plots for 1 model
#----------------------------------------------------------------------#
predPropCustom$PercentSATB2pos <- seq(0,100,10) # add col for

plotdf <- melt(predPropCustom[,c(1,2,5)], id="PercentSATB2pos")
colnames(plotdf)[2:3] <- c("CellType", "PredictedProportion")

# proportion plot
ggplot(plotdf, aes(x=PercentSATB2pos, y=PredictedProportion, colour=CellType))+
  geom_line()

# CETYGO box plot
ggplot(predPropCustom, aes(x=factor(0), y=CETYGO))+
         geom_boxplot()+
  geom_hline(yintercept = 0.07, colour="red", linetype="dashed")


# plots for all leave one out samples...
# take just the satb2pos proportions (and then again for neg), and plot against
# percent grouped by each sample



