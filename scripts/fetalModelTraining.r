#Emma Walker
#19/05/2023
#E.M.Walker@exeter.ac.uk

# create fetal model

# just for satb2+ and satb2- to begin with...

library(CETYGO)


setwd("/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/")

# load in fetal FANs data

load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/Fetal/FACS/2_normalised/normalisedBetas_FACS_Fetal.rdat")

SampleSheet$Cell_Type <- as.character(SampleSheet$Cell_Type)
SampleSheet$Cell_Type <- trimws(SampleSheet$Cell_Type)

# train model
# from github - https://github.com/ds420/CETYGO/blob/main/vignettes/QuantifyErrorInCellularCompositionEstimate.Rmd

# subset betas and pheno file to just satb2+ and -ve samples
table(SampleSheet$Cell_Type)

# change cell type names to remove special symbols
SampleSheet$Cell_Type <- gsub(" ", "", SampleSheet$Cell_Type, fixed = TRUE)
SampleSheet$Cell_Type <- gsub("SATB2+", "SATB2pos", SampleSheet$Cell_Type, fixed = T)
SampleSheet$Cell_Type <- gsub("SATB2-", "SATB2neg", SampleSheet$Cell_Type, fixed = T)


# subset to just satb2+ and satb- samples
SampleSheet <- SampleSheet[which(SampleSheet$Cell_Type %in% c("SATB2pos", "SATB2neg")),]
betas <- betas[,colnames(betas) %in% SampleSheet$Basename]

# make sure order of matches in pheno and betas
print(identical(as.character(SampleSheet$Basename), colnames(betas)))
#[1] TRUE



#Next, we select the sites to form the basis of the deconvolution.

customModel <- pickCompProbesMatrix(rawbetas = betas,
                                    cellTypes = unique(SampleSheet$Cell_Type),
                                    cellInd = SampleSheet$Cell_Type,
                                    numProbes = 100,
                                    probeSelect = "auto")


save(customModel, file="models/fetalBrain.rdat")

