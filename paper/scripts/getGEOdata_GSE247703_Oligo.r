##---------------------------------------------------------------------#
##
## Title: Get dc-hiOL data from GEO
##
## Methylation profiling of directly converted human fibroblasts from young, adult and old ## donors into oligodendrocytes (dc-hiOL) to study age-associated changes.
##
## Author: Emma Walker
##
## Date Created: 09/01/2024
## 
##---------------------------------------------------------------------#


#----------------------------------------------------------------------#
# USEFUL LINKS
#----------------------------------------------------------------------#

# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE247703
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8072064/#mmc1
# file:///C:/Users/EW367/Downloads/mmc1.pdf # supp info

#----------------------------------------------------------------------#
# ON COMMAND LINE
#----------------------------------------------------------------------#

# wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE247nnn/GSE247703/suppl//GSE247703_RAW.tar
# tar -xvf GSE247703_RAW.tar
# gunzip *.gz

library(GEOquery)

setwd("/lustre/projects/Research_Project-MRC190311/DNAm/cellLinePredictionCETGYO/data/")

GSE247703 <- getGEO('GSE247703',GSEMatrix=TRUE)

eData <-GSE247703[[1]]
df <- as.data.frame(pData(eData))

df$Basename <- substr(df$supplementary_file, nchar("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7899nnn/GSM7899124/suppl/")+1, nchar("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7899nnn/GSM7899124/suppl/GSM7899124_203168500073_R04C01") )

sampleSheet <- df %>% dplyr::select(Basename, source_name_ch1, "age_years:ch1", "gender:ch1")
colnames(sampleSheet) <- c("Basename", "Sample_ID", "Age", "Sex")

write.csv(sampleSheet, "GSE247703/QC/sampleSheet.csv")
