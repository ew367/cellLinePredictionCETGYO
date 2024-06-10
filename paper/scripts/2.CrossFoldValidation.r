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
##                                    14-20 weeks (most seperation for satb2+ satb2-)  
##        
##                   
##
##---------------------------------------------------------------------#

#----------------------------------------------------------------------#
# Notes
#----------------------------------------------------------------------#

# only include matched samples? - check how many are matched/unmatched


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
# Select samples to leave out in cross fold validation
#----------------------------------------------------------------------##

# n.b. some samples appear many times

# take 10 individuals that have one satb2+ and one satb2- sample

leaveOut <- c() 
for(i in unique(subPheno$Individual_ID)){
  indivCT <- sort(subPheno$Cell_Type[subPheno$Individual_ID == i], decreasing = F)
  if(length(indivCT) == 2){
    if(identical(c("SATB2neg","SATB2pos"), indivCT)){
      leaveOut <- c(leaveOut, i)
    }
  }
}

#randomly select 10 individuals
set.seed(1)
leaveOut <- sample(leaveOut,10)


#----------------------------------------------------------------------#
# Create model and run cross fold validation
#----------------------------------------------------------------------##

# for each in leaveOut
# create train and test data
# create model
# t
