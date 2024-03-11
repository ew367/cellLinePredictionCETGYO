##---------------------------------------------------------------------#
##
## Title: predict and correlate age of training data samples
##
## Purpose: predicted ages using fetal clock and cor with actual age
## 
##
## Author: Emma Walker
##
## Date Created: 08/03/2024
##
##---------------------------------------------------------------------##



#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

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


# source fetalClock
source("/lustre/projects/Research_Project-MRC190311/references/FetalClock/FetalClockFunction.R")


#----------------------------------------------------------------------#
# Run fetal clock
#----------------------------------------------------------------------#

agePred <- FetalClock(satb2Betas)
# [1] "Some probes of the coefficients are missing in the betas. We will need to impute values here - The final predictions will be less accurate"


class(agePred)

ageCorPlotdf <- as.data.frame(cbind(satb2$Age, agePred, satb2$Cell_Type))
colnames(ageCorPlotdf) <- c("Age", "predAge", "cellType")
ageCorPlotdf <- ageCorPlotdf %>% mutate(Age = as.numeric(as.character(Age)),
                                        predAge = as.numeric(as.character(predAge)))

# plot age vs pred age coloured by outlier status
ggplot(ageCorPlotdf, aes(x=Age, y=predAge, colour=cellType))+
  geom_point()

