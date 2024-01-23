
##---------------------------------------------------------------------#
##
## Title: plots of SCZ EWAS DMRs
##
## Purpose of script: plot the DMR results from the DMRcate package
##
## Author: Emma Walker
##
## Date Created: 24/11/23
##
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(ggplot2)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

projDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/SCZ"
dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/MRC"
resDir <- file.path(dataDir, "CETYGO/EWAS/DMRcate")
plotDir <- file.path(dataDir, "CETYGO/plots/")

cellType <- "NeuN+"
models <- c("modelLM", "nullCT")
mod=1
nCpGThreshold <- 3

manifestFile<-"/gpfs/mrc0/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv"
resFile <- file.path(resDir, paste0(cellType, models[mod], "DMRs.rdat"))
normData <- file.path(dataDir, "3_normalised/normalised.rdata")

plotFile <- paste0(plotDir, "DMRcate", cellType, models[mod],".pdf")


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(projDir)

load(resFile)
load(normData)

epicManifest<-fread(manifestFile, skip=7, fill=TRUE, data.table=F)
epicMan<-epicManifest[match(rownames(celltypeNormbeta), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]
QCmetrics<-QCmetrics[!is.na(QCmetrics$CCDNAmAge),]

# subset beta matrix to analysis samples
celltypeNormbeta<-celltypeNormbeta[,QCmetrics$Basename]
betas.anno <- cbind(celltypeNormbeta, as.data.frame(epicMan))
betas.anno$CHR <- paste0("chr", betas.anno$CHR)


#----------------------------------------------------------------------#
# PLOTS
#----------------------------------------------------------------------#

# doesn't work until Gviz updated
#groups <- c(Control="magenta", Schizophrenia="forestgreen")
#cols <- groups[as.character(QCmetrics$Phenotype)]

#pdf(paste0(plotDir, "DMRcateTest.pdf"), width=60, height=100)
#par(mfrow=c(1,1))
#DMR.plot(ranges=results.ranges, dmr=1, CpGs=celltypeNormbeta, phen.col=cols, what="Beta",
#arraytype="EPIC", genome="hg19", toscale=TRUE, plotmedians=TRUE)
#dev.off()

results.ranges <- results.ranges[which(results.ranges$no.cpgs >= nCpGThreshold),]

pdf(plotFile, height = 5, width = 15)

for(i in 1:length(results.ranges)){
  
  DMR <- as.data.frame(results.ranges)[i,]
  DMR$seqnames <- as.character(DMR$seqnames)

  # extract DMR cpg info 
  plotcpgs <- betas.anno[which(betas.anno$CHR == DMR$seqnames & betas.anno$MAPINFO >= DMR$start-1000 & betas.anno$MAPINFO <= DMR$end+1000),]
  row.names(plotcpgs) <- plotcpgs$MAPINFO
  
  DMRcpgs <- betas.anno[which(betas.anno$CHR == DMR$seqnames & betas.anno$MAPINFO >= DMR$start & betas.anno$MAPINFO <= DMR$end),]
  
  row.names(DMRcpgs) <- DMRcpgs$MAPINFO


  #check betas and pheno are in same order
  identical(colnames(DMRcpgs)[1:164], QCmetrics$Basename)

  plotdf <- cbind(t(plotcpgs[1:164]), QCmetrics %>% dplyr::select(Phenotype))
  plotdf <- reshape2::melt(plotdf)
  colnames(plotdf) <- c("Phenotype", "Location", "Methylation")
  plotdf$Location <- as.numeric(as.character(plotdf$Location))


  P <- ggplot(plotdf, aes(x=Location, y=Methylation, shape=Phenotype, color=Phenotype)) +
    geom_point()+
    geom_smooth()+
    annotate("rect", xmin=min(DMRcpgs$MAPINFO), xmax=max(DMRcpgs$MAPINFO), ymin=0, ymax=Inf, alpha=0.2)+
    ggtitle(paste0(DMR$seqnames, ":", DMR$start, "-", DMR$end, "  ", DMR$overlapping.genes))
  
  #print(P)


  ## just DMR

  #check betas and pheno are in same order
  identical(colnames(DMRcpgs)[1:164], QCmetrics$Basename)

  plotdf <- cbind(t(DMRcpgs[1:164]), QCmetrics %>% dplyr::select(Phenotype))
  plotdf <- reshape2::melt(plotdf)
  colnames(plotdf) <- c("Phenotype", "Location", "Methylation")
  plotdf$Location <- as.numeric(as.character(plotdf$Location))


  P2 <- ggplot(plotdf, aes(x=Location, y=Methylation, shape=Phenotype, color=Phenotype)) +
    geom_point()+
    geom_smooth(method = "lm")+
    annotate("rect", xmin=min(DMRcpgs$MAPINFO), xmax=max(DMRcpgs$MAPINFO), ymin=0, ymax=Inf, alpha=0.2)+
    ggtitle(paste0(DMR$seqnames, ":", DMR$start, "-", DMR$end, "  ", DMR$overlapping.genes))


  #print(P2)
  
  grid.arrange(P, P2, ncol=2)
  
}

dev.off()

