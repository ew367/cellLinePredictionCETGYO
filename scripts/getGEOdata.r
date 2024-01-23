

library(GEOquery)
# get GEO data GSE30338 (cell line)

load("data/GSE30338.rdat")

eData <- eList[[1]]
names(pData(eData))
experimentData(eData)


df <- as.data.frame(pData(eData))

dfAstro <- df[which(df$`cell type:ch1` == "astrocyte"),]

eList2 <- getGEOSuppFiles("GSE30338")

proc <- 
