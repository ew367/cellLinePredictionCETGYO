##---------------------------------------------------------------------#
##
## Title: Look at overlap between probes in the models
##
## Purpose of script: create venn diagrams for the probes included in the models
## 
## allFetalSatb2Model is all the Sat2b+ and Satb2- fetal samples available (79 samples)
## weeks20Satb2Model is restricted to just samples with age < 20 weeks (70 samples)
## 
##
## Author: Emma Walker
##
## Date Created: 26/01/2024
##
##---------------------------------------------------------------------##


#----------------------------------------------------------------------#
# 0. Define Parameters
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/paper/"


#----------------------------------------------------------------------#
# SET UP
#----------------------------------------------------------------------#

library(RColorBrewer)
library(VennDiagram)

setwd(projDir)

# load models
load("models/allFetalSatb2.rdat")
load("models/weeksTo20FetalSatb2.rdat")
load("models/weeks16to20FetalSatb2.rdat")


# load celltype EWAS res

res <- read.csv("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/EWAS/Cell/FACS_Cell_EWAS_sig_v2.csv")

res <- res[res$Sig.Fetal == T,]

#----------------------------------------------------------------------#
# venn diagram of model probes
#----------------------------------------------------------------------#

myCol <- brewer.pal(3, "Pastel2")

all <- rownames(allFetalSatb2model[[1]])
w16 <- rownames(weeks16to20Satb2Model[[1]])
w20 <- rownames(weeks20Satb2Model[[1]])


venn.diagram(
  x = list(all, w16, w20),
  category.names = c("All Samples" , "Weeks16to20" , "Weeks20"),
  filename = 'plots/modelProbes_venn_diagramm.png',
  fill = myCol,
  output=TRUE
)

