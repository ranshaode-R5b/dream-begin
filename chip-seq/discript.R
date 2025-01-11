setwd("/home/renshida/chipseq")
BiocManager::install(c("airway", "DESeq2", "edgeR", "limma", "ChIPpeakAnno", "ChIPseeker"))
BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")


#导入
rm(list = ls())
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
library(ChIPpeakAnno)
library(ChIPseeker)
#数据预处理
# 指定参考基因组
require(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene 
txdb

#读取 bed 文件
peak <- readPeakFile('BRD4_intersect_summits.bed')
keepChr <- !grepl("_",seqlevels(peak)) #去除带有“_”的染色体
seqlevels(peak,pruning.mode="coarse") <- seqlevels(peak)[keepChr]
#peaks 的基因组注释
peakAnno <- annotatePeak(peak, tssRegion = c(-3000, 3000), 
                         TxDb = txdb,annoDb = "org.Hs.eg.db")
#输出结果
write.table(peakAnno, 'peak.txt', sep = '\t', row.names = FALSE, quote = FALSE)
peakAnno_df <- as.data.frame(peakAnno)


#可视化
# 查看所有peaks在染色体的分布情况
covplot(peak,weightCol = "V5")
# 查看peaks在所有基因启动子附近的分布情况
promoter <- getPromoters(TxDb = txdb,upstream = 3000, downstream = 3000)
tagMatrix <- getTagMatrix(peak,windows = promoter)
tagHeatmap(tagMatrix)

# 注释可视化
plotAnnoBar(peakAnno)
vennpie(peakAnno)
plotAnnoPie(peakAnno)
plotDistToTSS(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)

# 查看peaks长度分布，只统计1000bp一下的peaks
peaksLength=abs(peakAnno_df$end-peakAnno_df$start)
peaksLength=peaksLength[peaksLength<1000]
hist(peaksLength, breaks = 5, col = "lightblue", xlim=c(0,10), xlab="peak length", main="histogram of peak length")

# 得到peaks相关基因之后再用去做富集分析即可
genes=unique(peakAnno_df$ENSEMBL)