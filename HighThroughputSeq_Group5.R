
### Load required R modules (if not available, reinstall as admin)
library(data.table)
library(Gviz)
library(GenomicFeatures)
library(org.Hs.eg.db)

### Set directory to work in
setwd("C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs")


# 09, IFIT2:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr10"  
myStart = 89302000 
myEnd = 89309000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049609-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049609-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049609.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr10",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049609_IFIT2.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 09, IFNB1:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr9"  
myStart = 21077100 
myEnd = 21078000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049609-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049609-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049609.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr9",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049609_IFNB1.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 09, ISG15
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr1"  
myStart = 1013400 
myEnd = 1014600 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049609-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049609-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049609.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr1",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049609_ISG15.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()


########################################################################################
########################################################################################
########################################################################################



# 10, IFIT2:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr10"  
myStart = 89302000 
myEnd = 89309000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049610-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049610-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049610.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr10",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049610_IFIT2.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 10, IFNB1:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr9"  
myStart = 21077100 
myEnd = 21078000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049610-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049610-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049610.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr9",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049610_IFNB1.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 10, ISG15
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr1"  
myStart = 1013400 
myEnd = 1014600 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049610-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049610-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049610.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr1",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049610_ISG15.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()


########################################################################################
########################################################################################
########################################################################################


# 11, IFIT2:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr10"  
myStart = 89302000 
myEnd = 89309000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049611-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049611-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049611.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr10",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049611_IFIT2.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 11, IFNB1:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr9"  
myStart = 21077100 
myEnd = 21078000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049611-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049611-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049611.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr9",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049611_IFNB1.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 11, ISG15
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr1"  
myStart = 1013400 
myEnd = 1014600 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049611-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049611-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049611.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr1",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049611_ISG15.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()


########################################################################################
########################################################################################
########################################################################################


# 12, IFIT2:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr10"  
myStart = 89302000 
myEnd = 89309000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049612-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049612-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049612.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr10",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049612_IFIT2.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 12, IFNB1:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr9"  
myStart = 21077100 
myEnd = 21078000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049612-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049612-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049612.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr9",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049612_IFNB1.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 12, ISG15
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr1"  
myStart = 1013400 
myEnd = 1014600 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049612-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049612-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049612.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr1",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049612_ISG15.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()


########################################################################################
########################################################################################
########################################################################################


# 15, IFIT2:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr10"  
myStart = 89302000 
myEnd = 89309000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049615-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049615-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049615.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr10",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049615_IFIT2.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 15, IFNB1:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr9"  
myStart = 21077100 
myEnd = 21078000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049615-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049615-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049615.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr9",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049615_IFNB1.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 15, ISG15
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr1"  
myStart = 1013400 
myEnd = 1014600 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049615-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049615-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049615.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr1",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049615_ISG15.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()


########################################################################################
########################################################################################
########################################################################################

# 16, IFIT2:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr10"  
myStart = 89302000 
myEnd = 89309000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049616-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049616-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049616.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr10",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049616_IFIT2.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 16, IFNB1:
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr9"  
myStart = 21077100 
myEnd = 21078000 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049616-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049616-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049616.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr9",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049616_IFNB1.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()



# 16, ISG15
### Specify genome name, chromosome number, and co-ordinates for mapping
myGenome = "hg38"
myChr = "chr1"  
myStart = 1013400 
myEnd = 1014600 


### Read in bedgraph files as simple text files (specify column headers)
bedFile1 <- fread('./SRR7049616-fwd.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile2 <- fread('./SRR7049616-rev.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))
bedFile3 <- fread('./SRR7049616.bedgraph', col.names = c('chromosome', 'start', 'end', 'value'))

### Determine the maximium depth value within the specified mapping co-ordinates
ChrData<-bedFile3[bedFile3$chromosome=="chr1",]
ChrData2<-ChrData[ChrData$start>myStart]
ChrData3<-ChrData2[ChrData2$end<myEnd,]
maxV<-max(ChrData3$value)

### Generate dataTracks - type ?DataTrack to see options
dataTrack1 <- DataTrack(range = bedFile1, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark blue", col = "black", background.title = "black", ylim=c(0,maxV))
dataTrack2 <- DataTrack(range = bedFile2, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark red", col = "black", background.title = "transparent", ylim=c(maxV,0)) 
dataTrack3 <- DataTrack(range = bedFile3, type = "a", chromosome=myChr, genome = myGenome, name = "Seq. Depth", fill = "dark green", col = "black", background.title = "transparent") 

### Generate genome and ideogram tracks
gtrack<-GenomeAxisTrack(col="black") 
itrack <- IdeogramTrack(genome = myGenome, chromosome = myChr)

### Read in UCSC genes and track 
ucscGenes <- UcscTrack(genome = myGenome, table="ncbiRefSeq", track = 'NCBI RefSeq', trackType="GeneRegionTrack",
                       chromosome=myChr, rstarts = "exonStarts", rends = "exonEnds",
                       gene = "name", symbol = 'name', transcript = "name",
                       strand = "strand", stacking = 'pack', showID = T, geneSymbol = T)

z <- ranges(ucscGenes)
mcols(z)$symbol <- mapIds(org.Hs.eg.db, gsub("\\.[1-9]$", "", mcols(z)$symbol), "SYMBOL","REFSEQ")
ucscGenes2 <- ucscGenes
ranges(ucscGenes2) <- z


### GENERATE PLOT ###
png(filename = "C:/Users/bmurp/OneDrive/Desktop/NYU/Classes/Fall_2020/Bioinformatics/02_NGS_Sequence_Alignment_and_Visualization/Assignments/Outputs/Figs/Gviz/SRR7049616_ISG15.png")
plotTracks(list(itrack,dataTrack1,dataTrack2,dataTrack3,gtrack,ucscGenes2), collapseTranscripts = "meta", transcriptAnnotation = "symbol", from = myStart, to = myEnd, sizes=c(0.05,0.2,0.2,0.2,0.1,0.2), type="hist", col.histogram=NA, cex.title=1, cex.axis=1, title.width=1.2)
dev.off()


########################################################################################
########################################################################################
########################################################################################

