
# check overlaps in DMRs and EWAS

dataDir <- "/lustre/projects/Research_Project-MRC190311/DNAm/MRC"
cellType <- "NeuN+"
manifestFile<-"/gpfs/mrc0/projects/Research_Project-MRC190311/references/EPICArray/MethylationEPIC_v-1-0_B5.csv"


EWASpath <- file.path(dataDir, "CETYGO/EWAS")
DMRpath <- file.path(EWASpath, "DMRcate")

EWASres <- file.path(EWASpath, paste0(cellType, "LM.rdata"))
DMRres <- file.path(DMRpath, paste0(cellType, "modelLMDMRs.rdat"))

load(EWASres)
load(DMRres)

# annotate res file
epicManifest<-fread(manifestFile, skip=7, fill=TRUE, data.table=F)
epicMan<-epicManifest[match(rownames(outtab), epicManifest$IlmnID),c("IlmnID","CHR","MAPINFO","UCSC_RefGene_Name","UCSC_RefGene_Group")]
outtab <- cbind(outtab, as.data.frame(epicMan))

outtab$CHR<-as.numeric(as.character(outtab$CHR))
outtab$MAPINFO<-as.numeric(as.character(outtab$MAPINFO))

# remove SNP probes
outtab<-outtab[-grep("rs", rownames(outtab)),]
outtab<-outtab[complete.cases(outtab),]

#change to numeric
outtab$CHR <- as.numeric(outtab$CHR)



# convert EWAS res to seq ranges object

outtab$location <- paste0("chr", outtab$CHR, ":", outtab$MAPINFO)
CpGs <- GRanges(outtab$location)

olap <- findOverlaps(CpGs, results.ranges)


