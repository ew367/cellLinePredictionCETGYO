#Emma Walker
#22/09/2023
#E.M.Walker@exeter.ac.uk

setwd("/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/")

# load predictor
load("models/fetalBrain.rdat")

fetalModelTest <- function(betas, project){
  
  ## identify which sites in the model are present in test data
  rInd<-rownames(betas)[rownames(betas) %in% rownames(customModel$coefEsts)]
  predPropCustom<-as.data.frame(projectCellTypeWithError(betas, customModel$coefEsts[rInd,]))
  
  
  # CETYGO score boxplot
  pCETYGO <- ggplot(predPropCustom, aes(factor(0), CETYGO))+
    geom_boxplot()+
    ggtitle(project)+
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank())
    
  
  # proportions boxplot
  plotdf <-as.data.frame(predPropCustom[,1:(ncol(predPropCustom)-2)])
  plotdf$Basename <- rownames(predPropCustom)
  plotdf <- reshape2::melt(plotdf)
  colnames(plotdf) <- c("Basename", "CellType", "Proportion")
  
  p <- ggplot(plotdf, aes(x=CellType, y = Proportion)) +
    geom_boxplot()+
    ggtitle(project)+
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
  

  
  returnList <- list(predPropCustom, pCETYGO, p)
  
}


################################### testing data


# DRI IPSC
load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/agedNeurons/2_normalised/Normalised_DRICambridge.rdat")

dri <- fetalModelTest(betas, "DRI")

pdf("plots/dri.pdf")
print(dri[[2]])
print(dri[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Sam THC - ShSy5

load("data/THC_rerun_normalised.rdat")

sam <- fetalModelTest(betas2, "THC-ShSy5")

pdf("plots/thc.pdf")
print(sam[[2]])
print(sam[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Jenny - neurons

load("data/DasenData_badprobesrm_Neuron.rdata")

jen <- fetalModelTest(data.3_N, "ipsc neurons")

pdf("plots/jennyNeurons.pdf")
print(jen[[2]])
print(jen[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# hypo - astrocytes

load("data/hypo_beta_pf_ds_hybremoved_bisremoved.RData")
  
astro <- fetalModelTest(betas, "hypo astrocytes")

pdf("plots/hypoAstro.pdf")
print(astro[[2]])
print(astro[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Jenny microglia

load("data/DasenData_badprobesrm_Microglia.rdata")

micro <- fetalModelTest(data.3_MG, "microglia")

pdf("plots/micro.pdf")
print(micro[[2]])
print(micro[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


#### Alice bulk fetal

betas <- read.csv("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/Fetal/Bulk/2_normalised/fetalBulk_EX3_23pcw_n91.csv")
row.names(betas) <- betas$X
betas <- betas[,-c(1)]

bulk <- fetalModelTest(betas, "bulk")
pdf("plots/fetalBulk.pdf")
print(bulk[[2]])
print(bulk[[3]])
dev.off()

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


#### Alice epic v2 sorted ####### need to look at this more closely!!!

load("/lustre/projects/Research_Project-MRC190311/DNAm/SFARI2023/3_normalised/SFARI2023_betas.rdat")

for(i in unique(SampleSheet$Cell_type)){
cellBetas <- betas[,colnames(betas) %in% SampleSheet$Basename[which(SampleSheet$Cell_type == i)]]

out <- fetalModelTest(cellBetas, paste0("epicV2 sorted ", i))
print(out[[2]])
print(out[[3]])
}

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# Leo data

# tolosa 

load("data/data.tolosa.Rdata")

tolosaIPSCbetas <- betas.tolosa[,colnames(betas.tolosa) %in% pheno.tolosa$Sample[which(pheno.tolosa$Cell_State == "iPSC")]]

tolosaIPSC <- fetalModelTest(tolosaIPSCbetas, "tolosa iPSC")
print(tolosaIPSC[[2]])
print(tolosaIPSC[[3]])


tolosaNbetas <- betas.tolosa[,colnames(betas.tolosa) %in% pheno.tolosa$Sample[which(pheno.tolosa$Cell_State == "Neuron")]]

tolosaN <- fetalModelTest(tolosaNbetas, "tolosa Neurons")
print(tolosaN[[2]])
print(tolosaN[[3]])

rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# kiselev 

kiselevipscbetas <- betas.kiselev[,colnames(betas.kiselev) %in% pheno.kiselev$full_name[which(pheno.kiselev$Cell_state == "iPSC")]]

kiselevipsc <- fetalModelTest(kiselevipscbetas, "kiselev iPSC")
print(kiselevipsc[[2]])
print(kiselevipsc[[3]])



kiselevNbetas <- betas.kiselev[,colnames(betas.kiselev) %in% pheno.kiselev$full_name[which(pheno.kiselev$Cell_state == "Neuron")]]

kiselevN <- fetalModelTest(kiselevNbetas, "kiselev Neurons")
print(kiselevN[[2]])
print(kiselevN[[3]])

kiselevNbetas <- betas.kiselev[,colnames(betas.kiselev) %in% pheno.kiselev$full_name[which(pheno.kiselev$Cell_state == "Neuron")]]

kiselevN <- fetalModelTest(kiselevNbetas, "kiselev Neurons")
print(kiselevN[[2]])
print(kiselevN[[3]])


rm(list=setdiff(ls(), c("customModel", "fetalModelTest")))


# kim 

##### hES ===

load("data/data.kim.Rdata")
pheno.kim$Cell_State <- as.character(pheno.kim$Cell_State)

kimHESbetas <- betas.kim[,colnames(betas.kim) %in% pheno.kim$Sample[which(pheno.kim$Cell_State == "hES")]]

kimhES <- fetalModelTest(kimHESbetas, "kim hES")
print(kimhES[[2]])
print(kimhES[[3]])


#Neuron

kimNbetas <- betas.kim[,colnames(betas.kim) %in% pheno.kim$Sample[which(pheno.kim$Cell_State == "Neuron")]]

kimN <- fetalModelTest(kimNbetas, "kim Neurons")
print(kimN[[2]])
print(kimN[[3]])


# NPC 


kimNPCbetas <- betas.kim[,colnames(betas.kim) %in% pheno.kim$Sample[which(pheno.kim$Cell_State == "NPC")]]

kimNPC <- fetalModelTest(kimNbetas, "kim NPC")
print(kimNPC[[2]])
print(kimNPC[[3]])


### TN samples from training datset

load("/gpfs/mrc0/projects/Research_Project-MRC190311/DNAm/Fetal/FACS/2_normalised/normalisedBetas_FACS_Fetal.rdat")

SampleSheet$Cell_Type <- as.character(SampleSheet$Cell_Type)
SampleSheet$Cell_Type <- trimws(SampleSheet$Cell_Type)


TNbetas <- betas[,colnames(betas) %in% SampleSheet$Basename[which(SampleSheet$Cell_Type == "TN")]]

TN <- fetalModelTest(TNbetas, "TN")
print(TN[[2]])
print(TN[[3]])


## Total Nuclei

TNbetas <- betas[,colnames(betas) %in% SampleSheet$Basename[which(SampleSheet$Cell_Type == "Total Nuclei")]]

TN <- fetalModelTest(TNbetas, "Total Nuclei")
print(TN[[2]])
print(TN[[3]])















