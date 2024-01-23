


## DMR plot function from Morteza

DMR.Plot <- function(DMR.chr,DMR.start,DMR.end,DMR.gene,CpG.Annotation,CpG.beta,Phenotype){
  DMR.chr=as.numeric(DMR.chr)
  DMR.start=as.numeric(DMR.start)
  DMR.end=as.numeric(DMR.end)
  CpG.Annotation$CHR = as.numeric(CpG.Annotation$CHR)
  CpG.Annotation=CpG.Annotation[which(CpG.Annotation$CHR==DMR.chr & CpG.Annotation$MAPINFO >= DMR.start & CpG.Annotation$MAPINFO <= DMR.end),]
  CpG.Annotation = cbind.data.frame(cpg=rownames(CpG.Annotation),chr=CpG.Annotation$CHR,start=CpG.Annotation$MAPINFO,end=CpG.Annotation$MAPINFO+8)
  CpG.Annotation = CpG.Annotation[order(CpG.Annotation$start,decreasing = F),]
  CpG.beta = CpG.beta[CpG.Annotation$cpg,]
  CpG.data = cbind.data.frame(CpG.Annotation,CpG.beta)
  
  cpg.range1<-makeGRangesFromDataFrame(CpG.data[,-1],    #convert into a Grange object
                                       keep.extra.columns=F,
                                       ignore.strand=T,
                                       seqinfo=NULL,
                                       seqnames.field=c("chr"),
                                       start.field=c("start"),
                                       end.field=c("start"),
                                       starts.in.df.are.0based=F)
  cpg.range2<-makeGRangesFromDataFrame(CpG.data[,-1],    #convert into a Grange object
                                       keep.extra.columns=T,
                                       ignore.strand=T,
                                       seqinfo=NULL,
                                       seqnames.field=c("chr"),
                                       start.field=c("start"),
                                       end.field=c("end"),
                                       starts.in.df.are.0based=F)
  atrack <- AnnotationTrack(cpg.range1, name = "CpG sites")
  
  dmr.range=GRanges(seqnames=paste0("chr",DMR.chr),ranges = IRanges(start = DMR.start,end = DMR.end,names = DMR.gene))
  dmrtrack <- AnnotationTrack(dmr.range, name = "CpG Island")
  
  itrack <- IdeogramTrack(genome = "hg19", chromosome = DMR.chr)
  gtrack <- GenomeAxisTrack()
  
  dTrack.ec 	<- DataTrack(cpg.range2, genome = "hg19", name = "Normalized beta value",groups = Phenotype, type = c("boxplot","a","g"), col=c("red", "blue"), legend = TRUE)
  
  plotTracks(list(itrack,gtrack,dmrtrack,atrack,dTrack.ec),from =DMR.start, to =DMR.end+8, background.title = "darkblue", main=DMR.gene)
}