



#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(doParallel)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

args<-commandArgs(trailingOnly = TRUE)
dataDir <- args[1]
#dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/MRC/"
cellType <- args[2]
#cellType <- "NeuN+"
mod <- args[3]
#mod <- "SexGroup"
#resDir <- args[4]
#resDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/MRC/CETYGO/EWAS/"


normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "CETYGO/EWAS/sexAnalysis/")
scriptsPath <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/SCZ/"

#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)

print(paste0("running EWAS on ", cellType, " cell type..."))
## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
identical(colnames(celltypeNormbeta),QCmetrics$Basename)

# # take top 100 rows for debugging
# celltypeNormbeta <- celltypeNormbeta[1:100,]

#----------------------------------------------------------------------#
# INTITATE PARALLEL ENV
#----------------------------------------------------------------------#
columns <- c('Estimate.','SE.','P.')
output.colnames <- paste0(columns, mod)
if(mod=="SexGroup"){ output.colnames <- c(paste0(columns,'Sex:Group'), paste0(columns,'Sex'), paste0(columns,'Group') ) }

print("Output colnames will be ...")
print(output.colnames)

nCores<-detectCores()
cl <- makeCluster(nCores-1)
registerDoParallel(cl)

print("Sourcing SCZ_EWAS_functions.r")
source(paste0(scriptsPath, "SCZ_EWAS_functions.r"))
print(paste("Model to be used:", mod))
ANOVA.model <- get(mod)



#5. Run EWAS =========================================================================================================
print(paste("Starting EWAS using", mod, "function"))
print(paste("Starting EWAS on", nrow(QCmetrics), "samples."))

clusterExport(cl, list(mod))
print("Starting EWAS")
res <- foreach(i=1:nrow(celltypeNormbeta), .combine=rbind, .verbose=F) %dopar%{
			ANOVA.model(row=as.numeric(celltypeNormbeta[i,]), QCmetrics)}
print("Finished EWAS.")

print("Finalising results table.")
rownames(res)<-rownames(celltypeNormbeta); colnames(res) <- output.colnames


#6. Save results =====================================================================================================
outfile <- paste0("SCZ_EWAS_", mod, "_", cellType, ".rds")
print(paste0("Saving results to ... ", resPath, outfile))
saveRDS(res, file=paste0(resPath, outfile))
