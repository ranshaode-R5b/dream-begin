setwd("/home/renshida/eqtm/data/human")
######illumina 450k hg19 -> hg38
library(rtracklayer)
library(GenomicRanges)


# 读取 CSV 注释文件（确保路径正确，并根据实际列名调整）
ann <- read.csv("humanmethylation450_15017482_v1-2.csv", skip = 7, stringsAsFactors = FALSE)
ann$MAPINFO <- as.numeric(ann$MAPINFO)
ann <- ann[!is.na(ann$MAPINFO), ]
ann$CHR <- ifelse(grepl("^chr", ann$CHR), ann$CHR, paste0("chr", ann$CHR))

# 构建 hg19 的 GRanges 对象（假设染色体列为 'chr'，位置列为 'pos'）
gr.hg19 <- GRanges(seqnames = ann$CHR, ranges = IRanges(start = ann$MAPINFO, width = 1))

names(gr.hg19) <- seq_len(length(gr.hg19))

# 导入链文件并执行 liftOver
chain <- import.chain("hg19ToHg38.over.chain")
gr.hg38_list <- liftOver(gr.hg19, chain)
# unlist 时保留名称（原始索引）
gr.hg38 <- unlist(gr.hg38_list)

# 新建两个列，初始化为 NA
ann$hg38_chr <- NA
ann$hg38_pos <- NA

# 利用转换结果的名称（即原始行号），将转换成功的结果填回原始数据
idx <- as.integer(names(gr.hg38))
ann$hg38_chr[idx] <- as.character(seqnames(gr.hg38))
ann$hg38_pos[idx] <- start(gr.hg38)

# 查看部分结果
head(ann[, c("CHR", "MAPINFO", "hg38_chr", "hg38_pos")])
ann_450k <- ann

#######EPIC v1.0
ann_EPIC <- read.csv("/home/renshida/eqtm/rawdata/human/infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7, stringsAsFactors = FALSE)




###450k  hg19
#install.packages("dplyr")  
library(readr)
library(dplyr)
library(AnnotationDbi)
# 读取文件

##########31282290
data1 <- read.table("31282290-cis-mir.txt", header = TRUE, sep = "")

# 仅保留第一列
#data1 <- data1 %>% select(1)
data1 <- merge(data1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data1)[17] <- "CpG_Gene"
colnames(data1)[18] <- "Feature"
colnames(data1)[19] <- "C_Chr"
colnames(data1)[20] <- "C_Pos"
colnames(data1)[6] <- "P-value"
data1$`P-value` <- as.numeric(data1$`P-value`)
data1$Beta <- as.numeric(data1$Beta)
data1$z <- abs(qnorm(data1$`P-value` / 2))
data1$SE <- abs(data1$Beta) / data1$z
data1 <- data1[,-21]
write.table(data1, file = "/home/renshida/eqtm/data/data_cpg/31282290-cis-mir.txt", sep = "\t", row.names = FALSE)




library(ggplot2)
feature_counts <- table(data1$Feature)
feature_df <- as.data.frame(feature_counts)
colnames(feature_df) <- c("Feature", "Count")
ggplot(feature_df, aes(x = Feature, y = Count)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_minimal() +
  labs(title = "基因组特征分布",
       x = "基因组特征",
       y = "计数") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))






########37233989
data3 <- read.table("37233989-cis-vc.txt", header = TRUE, sep = "")

# 仅保留第一列
#data3 <- data3 %>% select(1)
data3 <- merge(data3, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data3)[17] <- "CpG_Gene"
colnames(data3)[18] <- "Feature"
colnames(data3)[19] <- "C_Chr"
colnames(data3)[20] <- "C_Pos"
colnames(data3)[6] <- "P-value"
data3$`P-value` <- as.numeric(data3$`P-value`)
data3$Beta <- as.numeric(data3$Beta)
data3$z <- abs(qnorm(data3$`P-value` / 2))
data3$SE <- abs(data3$Beta) / data3$z
data3 <- data3[,-21]
write.table(data3, file = "/home/renshida/eqtm/data/data_cpg/37233989-cis-vc.txt", sep = "\t", row.names = FALSE)


data3.1 <- read.table("37233989-cis-ve.txt", header = TRUE, sep = "")
# 仅保留第一列
#data3.1 <- data3.1 %>% select(1)
data3.1 <- merge(data3.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data3.1)[17] <- "CpG_Gene"
colnames(data3.1)[18] <- "Feature"
colnames(data3.1)[19] <- "C_Chr"
colnames(data3.1)[20] <- "C_Pos"
colnames(data3.1)[6] <- "P-value"
data3.1$`P-value` <- as.numeric(data3.1$`P-value`)
data3.1$Beta <- as.numeric(data3.1$Beta)
data3.1$z <- abs(qnorm(data3.1$`P-value` / 2))
data3.1$SE <- abs(data3.1$Beta) / data3.1$z
data3.1 <- data3.1[,-21]
write.table(data3.1, file = "/home/renshida/eqtm/data/data_cpg/37233989-cis-ve.txt", sep = "\t", row.names = FALSE)




#########39393618
data4 <- read.table("39393618-cis.txt",header = TRUE,sep="")
#data4 <- data4 %>% select(1)
data4 <- merge(data4, 
                 ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data4)[17] <- "CpG_Gene"
colnames(data4)[18] <- "Feature"
colnames(data4)[19] <- "C_Chr"
colnames(data4)[20] <- "C_Pos"
colnames(data4)[6] <- "P-value"
data4$`P-value` <- as.numeric(data4$`P-value`)
data4$Beta <- as.numeric(data4$Beta)
data4$z <- abs(qnorm(data4$`P-value` / 2))
data4$SE <- abs(data4$Beta) / data4$z
data4 <- data4[,-21]
write.table(data4, file = "/home/renshida/eqtm/data/data_cpg/39393618-cis.txt", sep = "\t", row.names = FALSE)
















