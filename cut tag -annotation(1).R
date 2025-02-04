library(devtools)
library(TxDb.Hsapiens.UCSC.hg19.knownGene) 
library(ChIPpeakAnno) 
library(ChIPseeker)
require(ChIPseeker)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(clusterProfiler) 
RUNX1_bedPeaksFile   = "RUNX1_peaks.xls"; 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
RUNX1_peak <- readPeakFile(RUNX1_bedPeaksFile)  
keepChr= !grepl('_',seqlevels(RUNX1_peak))
seqlevels(RUNX1_peak,pruning.mode="coarse") <- seqlevels(RUNX1_peak)[keepChr]
RUNX1_peakAnno <- annotatePeak(RUNX1_peak, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db") 
TP63_bedPeaksFile   = "TP63_peaks.xls"; 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
TP63_peak <- readPeakFile(TP63_bedPeaksFile)  
keepChr= !grepl('_',seqlevels(TP63_peak))
seqlevels(TP63_peak,pruning.mode="coarse") <- seqlevels(TP63_peak)[keepChr]
TP63_peakAnno <- annotatePeak(TP63_peak, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db") 
GLI2_bedPeaksFile   = "GLI2_peaks.xls"; 
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
GLI2_peak <- readPeakFile(GLI2_bedPeaksFile)  
keepChr= !grepl('_',seqlevels(GLI2_peak))
seqlevels(GLI2_peak,pruning.mode="coarse") <- seqlevels(GLI2_peak)[keepChr]
GLI2_peakAnno <- annotatePeak(GLI2_peak, tssRegion=c(-3000, 3000), 
                         TxDb=txdb, annoDb="org.Hs.eg.db") 