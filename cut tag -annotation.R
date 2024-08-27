library(devtools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
library(ChIPpeakAnno) 
library(ChIPseeker)
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(clusterProfiler) 
bedPeaksFile   = "RUNX1_peaks.xls"; 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peak <- readPeakFile(bedPeaksFile)  
keepChr= !grepl('_',seqlevels(peak))
seqlevels(peak,pruning.mode="coarse") <- seqlevels(peak)[keepChr]
peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db") 
plotAnnoPie(peakAnno)
