##---------------------------------------------------------------------#
##
## Title: Look at overlap between probes in the models
##
## Purpose of script: create venn diagrams for the probes included in the models
## 
## allFetalSatb2Modelx is all the Sat2b+ and Satb2- fetal samples available 
## (79 samples)
## weeks20Satb2Model is restricted to just samples with age < 20 weeks 
## (70 samples)
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
load("models/allFetalSatb2Any.rdat")
load("models/allFetalSatb2Auto.rdat")
load("models/weeksTo20FetalSatb2.rdat")
load("models/weeks16to20FetalSatb2.rdat")


# load celltype EWAS res
res <- read.csv("/lustre/projects/Research_Project-MRC190311/DNAm/Lifecourse1/FANS/3_analysis/EWAS/FACS_Cell_EWAS_fetal_adult_anno.csv")
res <- res[res$Sig.Fetal == T,]

# extract probeIDs

auto <- rownames(allFetalSatb2modelAuto[[1]])
any <- rownames(allFetalSatb2modelAny[[1]])
w16 <- rownames(weeks16to20Satb2Model[[1]])
w20 <- rownames(weeks20Satb2Model[[1]])
resID <- res$IlmnID

#----------------------------------------------------------------------#
# venn diagrams of model probes
#----------------------------------------------------------------------#

#all models

myCol <- brewer.pal(4, "Pastel2")

venn.diagram(
  x = list(auto, any, w16, w20),
  category.names = c("All Samples auto", "All Samples any", "Weeks16to20" , "Weeks20"),
  filename = 'plots/vennDiagrams/modelProbes_venn_diagramm.png',
  fill = myCol,
  output=TRUE
)


# res/auto/any

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(auto, any, res),
  category.names = c("All Samples auto", "All Samples any", "EWAS"),
  filename = 'plots/vennDiagrams/ResAutoAnymodelProbes_venn_diagramm.png',
  fill = myCol,
  output=TRUE
)



# res top 200/auto/any

res200ID <- res$IlmnID[order(res$Fetal.P.Cell)[1:200]]

myCol <- brewer.pal(3, "Pastel2")

venn.diagram(
  x = list(auto, any, res200ID),
  category.names = c("All Samples auto", "All Samples any", "EWAS"),
  filename = 'plots/vennDiagrams/Res200AutoAnymodelProbes_venn_diagramm.png',
  fill = myCol,
  output=TRUE
)





