##---------------------------------------------------------------------#
##
## Title: DMRFF analysis
##
## Purpose of script: perform regional analysis of 
## schizophrenia vs controls EWAS results
## Author: Emma Walker
##
## Date Created: 07-11-2023
##
##---------------------------------------------------------------------#

######## N.B. THIS SCRIPT IS NOT RUNNING AS A JOB SUBMISSON

#----------------------------------------------------------------------#
# LOAD PACKAGES
#----------------------------------------------------------------------#

library(data.table)
library(dmrff)

#----------------------------------------------------------------------#
# DEFINE PARAMETERS
#----------------------------------------------------------------------#

#args<-commandArgs(trailingOnly = TRUE)
#dataDir <- args[1]
dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/MRC"
#cellType <- args[2]
cellType <- "Sox10+"

maxGap <- 500

normData<-file.path(dataDir, "3_normalised/normalised.rdata")
resPath<-file.path(dataDir, "CETYGO/EWAS/")
resFile<-file.path(resPath, paste0(cellType, "LM.rdata"))
manifestFile<-"/gpfs/mrc0/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv"


#----------------------------------------------------------------------#
# LOAD AND PREPARE DATA
#----------------------------------------------------------------------#

setwd(dataDir)
load(normData)
load(resFile)

# annotate results file
epicManifest<-fread(manifestFile, skip=7, fill=TRUE, data.table=F)
epicMan<-epicManifest[match(rownames(outtab), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
outtab <- cbind(outtab, as.data.frame(epicMan))

# remove sex chrs
outtab$CHR<-as.numeric(as.character(outtab$CHR))
outtab$MAPINFO<-as.numeric(as.character(outtab$MAPINFO))

# remove SNP probes
outtab<-outtab[-grep("rs", rownames(outtab)),]
outtab<-outtab[complete.cases(outtab),]

#change to numeric
outtab$CHR <- as.numeric(outtab$CHR)

## subset to cell type samples
QCmetrics<-QCmetrics[which(QCmetrics$Cell.type == cellType),]

# filter to schizophrenia and control only
QCmetrics<-QCmetrics[QCmetrics$Phenotype %in% c("Schizophrenia", "Control"),]

# subset beta matrix to analysis samples and cpgs in outtab
celltypeNormbeta<-celltypeNormbeta[rownames(celltypeNormbeta) %in% rownames(outtab),QCmetrics$Basename]


#----------------------------------------------------------------------#
# RUN DMRFF
#----------------------------------------------------------------------#

## DMRFF requirements

# data frame stats which has the following columns:
#estimate (regression coefficient),
#se (standard error of the coefficient),
#p.value,
#chr (chromosome of the CpG site),
#pos (position of the CpG site on the chromosome).

#DNA methylation dataset in R as matrix methylation for which rows correspond to CpG sites and columns to samples. The DNA methylation levels are necessary for dmrff to calculate and adjust for dependencies between CpG sites.

print(paste0("Running DMRFF on ", cellType, " modelLM results..."))

dmrsLM <- dmrff(estimate=outtab$SCZ_coeff,
              se=outtab$SCZ_SE,
              p.value=outtab$SCZ_P,
              methylation=celltypeNormbeta,
              chr=outtab$CHR,
              pos=outtab$MAPINFO,
              maxgap=maxGap,
              verbose=T)

save(dmrsLM, file = file.path(paste0(resPath, cellType,"LMdmrs.rdata")))
#save(dmrsLM, file = file.path(paste0(resPath, cellType,"nullCTdmrs.rdata")))

if (cellType == "Double-" | cellType == "NeuN+"){
  
  print(paste0("Running DMRFF on ", cellType, " nullCT results..."))
  
  dmrsNullCT <- dmrff(estimate=outtab$nullCT_SCZ_coeff,
                  se=outtab$nullCT_SCZ_SE,
                  p.value=outtab$nullCT_SCZ_P,
                  methylation=celltypeNormbeta,
                  chr=outtab$CHR,
                  pos=outtab$MAPINFO,
                  maxgap=maxGap,
                  verbose=T)
  
  save(dmrsNullCT, file = file.path(paste0(resPath, cellType,"nullCTdmrs.rdata")))

} else if (cellType == "Sox10+"){
  
  print("No nullCT model for Sox10+ data")
  
} else {
  
  stop("Cell Type not recognised....")
  
}

print("job finished")


#----------------------------------------------------------------------#
# ANALYSE OUTPUT
#----------------------------------------------------------------------#

dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/MRC"
cellType <- "Double-"
resPath<-file.path(dataDir, "CETYGO/EWAS/")

load(file.path(paste0(resPath, cellType,"LMdmrs.rdata")))

# NOTE for SOX10 the object is called dmrsLM but is actually cellCTdmrs as no cell proportion data for this fraction

# from Emma's script

dmrsLM <- dmrsLM[which(dmrsLM$p.adjust < 0.05 & dmrsLM$n > 2),] 
dim(dmrsLM) ## 1 for double-, 
dmrsLM

## get the CpG sites which are in each region 
sites <- dmrff.sites(dmrsLM, chr=outtab$CHR, outtab$MAPINFO)
sites$location <- paste0(sites$chr, "_", sites$pos)
outtab$location <- paste0(outtab$CHR, "_", outtab$MAPINFO)


sites <- dplyr::left_join(sites, outtab %>% dplyr::select("UCSC_RefGene_Name", "SCZ_coeff", "SCZ_P", "location"))
sites <- dplyr::left_join(sites, dmrsLM %>% dplyr::select("start","end","z","p.adjust", "chr"))

outFile <- paste0(resPath, cellType, "DMR.csv")
write.csv(dmrsLM, outFile )
outFile <- paste0(resPath, cellType, "DMRsites.csv")
write.csv(sites,outFile)

## to look at sites with the max regions

#max.region <- which.max(dmrs$n)
#kable(sites[which(sites$region == max.region),],row.names=F)

dmrff.sites <- function(regions, chr, pos) {
  stopifnot(is.data.frame(regions) && all(c("chr","start","end") %in% colnames(regions)))
  stopifnot(is.vector(chr))
  stopifnot(is.vector(pos))
  stopifnot(length(chr) == length(pos))
  
  sites <- data.frame(site=1:length(chr),
                      chr=chr,
                      pos=pos)
  
  members <- region.members(regions, sites)
  
  idx <- unlist(members)
  if (length(idx) > 0)
    data.frame(region=rep(1:nrow(regions), sapply(members, length)),
               sites[idx,],
               stringsAsFactors=F)
  else
    NULL
}


region.members <- function(intervals, positions) {
  stopifnot(is.data.frame(intervals))
  stopifnot(all(c("chr","start","end") %in% colnames(intervals)))
  stopifnot(is.data.frame(positions))
  stopifnot(all(c("chr","pos") %in% colnames(positions)))
  
  events <- rbind(data.frame(chr=as.character(intervals$chr),
                             pos=intervals$start,
                             type="start",
                             id=1:nrow(intervals)),
                  data.frame(chr=as.character(intervals$chr),
                             pos=intervals$end,
                             type="end",
                             id=1:nrow(intervals)),
                  data.frame(chr=as.character(positions$chr),
                             pos=positions$pos,
                             type="pos",
                             id=1:nrow(positions)))
  events <- events[order(events$chr, events$pos, sign(events$type!="start"), sign(events$type=="end"), decreasing=F),]
  start.idx <- which(events$type == "start")
  end.idx <- which(events$type == "end")
  intervals$event.start.idx <- start.idx[match(1:nrow(intervals), events$id[start.idx])] + 1
  intervals$event.end.idx <- end.idx[match(1:nrow(intervals), events$id[end.idx])] - 1
  events$id[which(events$type != "pos")] <- NA    
  lapply(1:nrow(intervals), function(i) {
    if (intervals$event.start.idx[i] >= intervals$event.end.idx[i]) integer(0)
    na.omit(events$id[intervals$event.start.idx[i]:intervals$event.end.idx[i]])
  })         
}



dmrff.manhattan.plot(dmrsLM, chr=outtab$CHR, outtab$MAPINFO, title = "Manhattan plot")




