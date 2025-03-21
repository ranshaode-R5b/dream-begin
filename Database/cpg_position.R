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




#######37271320
data2 <- read.table("37271320-cis.txt", header = TRUE, sep = "")

# 仅保留第一列
#data2 <- data2 %>% select(1)
data2 <- merge(data2, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2)[17] <- "CpG_Gene"
colnames(data2)[18] <- "Feature"
colnames(data2)[19] <- "C_Chr"
colnames(data2)[20] <- "C_Pos"
colnames(data2)[6] <- "P-value"
data2$`P-value` <- as.numeric(data2$`P-value`)
data2$Beta <- as.numeric(data2$Beta)
data2$z <- abs(qnorm(data2$`P-value` / 2))
data2$SE <- abs(data2$Beta) / data2$z
data2 <- data2[,-21]
data2[1294:1295, 19] <- "chr22"
data2[1294:1295, 20] <- 24372133
data2[3712:3714, 19] <- "chr22"
data2[3712:3714, 20] <- 24372951
data2[6024:6025,19] <- "chr15"
data2[6024:6025,20] <- 82959442
data2[6261,19] <- "chr17"
data2[6261,20] <- 34611341
data2[6421:6426,19] <- "chr22"
data2[6421:6426,20] <- 24372958
write.table(data2, file = "/home/renshida/eqtm/data/data_cpg/37271320-cis.txt", sep = "\t", row.names = FALSE)


data2.1 <- read.table("37271320-trans1.txt", header = TRUE, sep = "")
# 仅保留第一列
#data2.1 <- data2.1 %>% select(1)
data2.1 <- merge(data2.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2.1)[17] <- "CpG_Gene"
colnames(data2.1)[18] <- "Feature"
colnames(data2.1)[19] <- "C_Chr"
colnames(data2.1)[20] <- "C_Pos"
colnames(data2.1)[6] <- "P-value"
data2.1$`P-value` <- as.numeric(data2.1$`P-value`)
data2.1$Beta <- as.numeric(data2.1$Beta)
data2.1$z <- abs(qnorm(data2.1$`P-value` / 2))
data2.1$SE <- abs(data2.1$Beta) / data2.1$z
data2.1 <- data2.1[,-21]
write.table(data2.1, file = "/home/renshida/eqtm/data/data_cpg/37271320-trans1.txt", sep = "\t", row.names = FALSE)

data2.2 <- read.table("37271320-trans2.txt", header = TRUE, sep = "")
# 仅保留第一列
#data2.2 <- data2.2 %>% select(1)
data2.2 <- merge(data2.2, 
                 ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2.2)[17] <- "CpG_Gene"
colnames(data2.2)[18] <- "Feature"
colnames(data2.2)[19] <- "C_Chr"
colnames(data2.2)[20] <- "C_Pos"
colnames(data2.2)[6] <- "P-value"
data2.2$`P-value` <- as.numeric(data2.2$`P-value`)
data2.2$Beta <- as.numeric(data2.2$Beta)
data2.2$z <- abs(qnorm(data2.2$`P-value` / 2))
data2.2$SE <- abs(data2.2$Beta) / data2.2$z
data2.2 <- data2.2[,-21]
write.table(data2.2, file = "/home/renshida/eqtm/data/data_cpg/37271320-trans2.txt", sep = "\t", row.names = FALSE)


data2.3 <- read.table("37271320-trans3.txt", header = TRUE, sep = "")
# 仅保留第一列
#data2.3 <- data2.3 %>% select(1)
data2.3 <- merge(data2.3, 
                 ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2.3)[17] <- "CpG_Gene"
colnames(data2.3)[18] <- "Feature"
colnames(data2.3)[19] <- "C_Chr"
colnames(data2.3)[20] <- "C_Pos"
colnames(data2.3)[6] <- "P-value"
data2.3$`P-value` <- as.numeric(data2.3$`P-value`)
data2.3$Beta <- as.numeric(data2.3$Beta)
data2.3$z <- abs(qnorm(data2.3$`P-value` / 2))
data2.3$SE <- abs(data2.3$Beta) / data2.3$z
data2.3 <- data2.3[,-21]
write.table(data2.3, file = "/home/renshida/eqtm/data/data_cpg/37271320-trans3.txt", sep = "\t", row.names = FALSE)


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



#######31076557
data5 <- read.table("31076557-cis.txt",header = TRUE,sep="")
#data5 <- data5 %>% select(1)
data5 <- merge(data5, 
               ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data5)[17] <- "CpG_Gene"
colnames(data5)[18] <- "Feature"
colnames(data5)[19] <- "C_Chr"
colnames(data5)[20] <- "C_Pos"
colnames(data5)[6] <- "P-value"
data5$`P-value` <- as.numeric(data5$`P-value`)
data5$Beta <- as.numeric(data5$Beta)
data5$z <- abs(qnorm(data5$`P-value` / 2))
data5$SE <- abs(data5$Beta) / data5$z
data5 <- data5[,-21]
write.table(data5, file = "/home/renshida/eqtm/data/data_cpg/31076557-cis.txt", sep = "\t", row.names = FALSE)


########35710981
data6 <- read.table("35710981-cis.txt",header = TRUE,sep="")
#data6 <- data6 %>% select(1)
data6 <- merge(data6, 
               ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data6)[17] <- "CpG_Gene"
colnames(data6)[18] <- "Feature"
colnames(data6)[19] <- "C_Chr"
colnames(data6)[20] <- "C_Pos"
colnames(data6)[6] <- "P-value"
data6$`P-value` <- as.numeric(data6$`P-value`)
data6$Beta <- as.numeric(data6$Beta)
data6$z <- abs(qnorm(data6$`P-value` / 2))
data6$SE <- abs(data6$Beta) / data6$z
data6 <- data6[,-21]
write.table(data6, file = "/home/renshida/eqtm/data/data_cpg/35710981-cis.txt", sep = "\t", row.names = FALSE)

#######31791327
data7 <- read.table("31791327-cis-LL.txt",header = TRUE,sep="")

#data7 <- data7 %>% select(1)
data7 <- merge(data7, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7)[17] <- "CpG_Gene"
colnames(data7)[18] <- "Feature"
colnames(data7)[19] <- "C_Chr"
colnames(data7)[20] <- "C_Pos"
colnames(data7)[6] <- "P-value"
write.table(data7, file = "/home/renshida/eqtm/data/data_cpg/31791327-cis-LL.txt", sep = "\t", row.names = FALSE)

data7.1 <- read.table("31791327-cis-LLS.txt",header = TRUE,sep="")
#data7.1 <- data7.1 %>% select(1)
data7.1 <- merge(data7.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7.1)[17] <- "CpG_Gene"
colnames(data7.1)[18] <- "Feature"
colnames(data7.1)[19] <- "C_Chr"
colnames(data7.1)[20] <- "C_Pos"
colnames(data7.1)[6] <- "P-value"
write.table(data7.1, file = "/home/renshida/eqtm/data/data_cpg/31791327-cis-LLS.txt", sep = "\t", row.names = FALSE)


data7.2 <- read.table("31791327-cis-NTR.txt",header = TRUE,sep="")
#data7.2 <- data7.2 %>% select(1)
data7.2 <- merge(data7.2, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7.2)[17] <- "CpG_Gene"
colnames(data7.2)[18] <- "Feature"
colnames(data7.2)[19] <- "C_Chr"
colnames(data7.2)[20] <- "C_Pos"
colnames(data7.2)[6] <- "P-value"
write.table(data7.2, file = "/home/renshida/eqtm/data/data_cpg/31791327-cis-NTR.txt", sep = "\t", row.names = FALSE)



data7.3 <- read.table("31791327-cis-RS.txt",header = TRUE,sep="")
#data7.3 <- data7.3 %>% select(1)
data7.3 <- merge(data7.3, 
                 ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7.3)[17] <- "CpG_Gene"
colnames(data7.3)[18] <- "Feature"
colnames(data7.3)[19] <- "C_Chr"
colnames(data7.3)[20] <- "C_Pos"
colnames(data7.3)[6] <- "P-value"
write.table(data7.3, file = "/home/renshida/eqtm/data/data_cpg/31791327-cis-RS.txt", sep = "\t", row.names = FALSE)

###############35302492
# data8 <- read.table("35302492-cis-adj.txt",header = TRUE,sep="")
# 
# data8 <- data8 %>% select(1)
# data8 <- merge(data8, 
#                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
#                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
# colnames(data8)[2] <- "CpG_Gene"
# colnames(data8)[3] <- "Feature"
# colnames(data8)[4] <- "Chr"
# colnames(data8)[5] <- "Pos"


data8.1 <- read.table("35302492-cis-unadj.txt",header = TRUE,sep="")

#data8.1 <- data8.1 %>% select(1)
data8.1 <- merge(data8.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data8.1)[17] <- "CpG_Gene"
colnames(data8.1)[18] <- "Feature"
colnames(data8.1)[19] <- "C_Chr"
colnames(data8.1)[20] <- "C_Pos"
write.table(data8.1, file = "/home/renshida/eqtm/data/data_cpg/35302492-cis-unadj.txt", sep = "\t", row.names = FALSE)



######37598461
data9 <- read.table("37598461-cis.txt",header = TRUE,sep="")
#data9 <- data9 %>% select(1)
data9 <- merge(data9, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data9)[17] <- "CpG_Gene"
colnames(data9)[18] <- "Feature"
colnames(data9)[19] <- "C_Chr"
colnames(data9)[20] <- "C_Pos"
colnames(data9)[6] <- "P-value"
data9$`P-value` <- as.numeric(data9$`P-value`)
data9$Beta <- as.numeric(data9$Beta)
data9$z <- abs(qnorm(data9$`P-value` / 2))
data9$SE <- abs(data9$Beta) / data9$z
data9 <- data9[,-21]
write.table(data9, file = "/home/renshida/eqtm/data/data_cpg/37598461-cis.txt", sep = "\t", row.names = FALSE)



#######33092652
data10 <- read.table("33092652-cis-trans-kora.txt",header = TRUE,sep="")
#colnames(data10)[1] <- "CpG"
#data10 <- data10 %>% select(1)
data10 <- merge(data10, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data10)[17] <- "CpG_Gene"
colnames(data10)[18] <- "Feature"
colnames(data10)[19] <- "C_Chr"
colnames(data10)[20] <- "C_Pos"
colnames(data10)[6] <- "P-value"
data10$`P-value` <- as.numeric(data10$`P-value`)
data10$Beta <- as.numeric(data10$Beta)
data10$z <- abs(qnorm(data10$`P-value` / 2))
data10$SE <- abs(data10$Beta) / data10$z
data10 <- data10[,-21]
write.table(data10, file = "/home/renshida/eqtm/data/data_cpg/33092652-cis-trans-kora.txt", sep = "\t", row.names = FALSE)

data10.1 <- read.table("33092652-cis-trans-rot.txt",header = TRUE,sep="")
#colnames(data10.1)[1] <- "CpG"
#data10.1 <- data10.1 %>% select(1)
data10.1 <- merge(data10.1, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data10.1)[17] <- "CpG_Gene"
colnames(data10.1)[18] <- "Feature"
colnames(data10.1)[19] <- "C_Chr"
colnames(data10.1)[20] <- "C_Pos"
colnames(data10.1)[6] <- "P-value"
data10.1$`P-value` <- as.numeric(data10.1$`P-value`)
data10.1$Beta <- as.numeric(data10.1$Beta)
data10.1$z <- abs(qnorm(data10.1$`P-value` / 2))
data10.1$SE <- abs(data10.1$Beta) / data10.1$z
data10.1 <- data10.1[,-21]
write.table(data10.1, file = "/home/renshida/eqtm/data/data_cpg/33092652-cis-trans-rot.txt", sep = "\t", row.names = FALSE)



#######29914364
data11 <- read.table("29914364-GTP.txt",header = TRUE,sep="")
#data11 <- data11 %>% select(1)
data11 <- merge(data11, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data11)[17] <- "CpG_Gene"
colnames(data11)[18] <- "Feature"
colnames(data11)[19] <- "C_Chr"
colnames(data11)[20] <- "C_Pos"
colnames(data11)[6] <- "P-value"
data11$`P-value` <- as.numeric(data11$`P-value`)
data11$Beta <- as.numeric(data11$Beta)
data11$z <- abs(qnorm(data11$`P-value` / 2))
data11$SE <- abs(data11$Beta) / data11$z
data11 <- data11[,-21]
write.table(data11, file = "/home/renshida/eqtm/data/data_cpg/29914364-GTP.txt", sep = "\t", row.names = FALSE)

data11.1 <- read.table("29914364-MESA.txt",header = TRUE,sep="")
#data11.1 <- data11.1 %>% select(1)
data11.1 <- merge(data11.1, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data11.1)[17] <- "CpG_Gene"
colnames(data11.1)[18] <- "Feature"
colnames(data11.1)[19] <- "C_Chr"
colnames(data11.1)[20] <- "C_Pos"
colnames(data11.1)[6] <- "P-value"
data11.1$`P-value` <- as.numeric(data11.1$`P-value`)
data11.1$Beta <- as.numeric(data11.1$Beta)
data11.1$z <- abs(qnorm(data11.1$`P-value` / 2))
data11.1$SE <- abs(data11.1$Beta) / data11.1$z
data11.1 <- data11.1[,-21]
write.table(data11.1, file = "/home/renshida/eqtm/data/data_cpg/29914364-MESA.txt", sep = "\t", row.names = FALSE)



#############35426765
data12 <- read.table("35426765-cis.txt",header = TRUE,sep = "")
#colnames(data12)[1] <- "CpG"
#data12 <- data12 %>% select(1)
data12 <- merge(data12, 
               ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data12)[17] <- "CpG_Gene"
colnames(data12)[18] <- "Feature"
colnames(data12)[19] <- "C_Chr"
colnames(data12)[20] <- "C_Pos"
colnames(data12)[6] <- "P-value"
write.table(data12, file = "/home/renshida/eqtm/data/data_cpg/35426765-cis.txt", sep = "\t", row.names = FALSE)


###########38438351
data13 <- read.table("38438351-cis.txt",header = TRUE,sep = "")
#colnames(data13)[1] <- "CpG"
#data13 <- data13 %>% select(1)
data13 <- merge(data13, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data13)[17] <- "CpG_Gene"
colnames(data13)[18] <- "Feature"
colnames(data13)[19] <- "C_Chr"
colnames(data13)[20] <- "C_Pos"
colnames(data13)[6] <- "P-value"
write.table(data13, file = "/home/renshida/eqtm/data/data_cpg/38438351-cis.txt", sep = "\t", row.names = FALSE)



######32493484
data14 <- read.table("32493484-cis.txt",header = TRUE,sep="")
#colnames(data14)[1] <- "CpG"
#data14 <- data14 %>% select(1)
data14 <- merge(data14, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data14)[17] <- "CpG_Gene"
colnames(data14)[18] <- "Feature"
colnames(data14)[19] <- "C_Chr"
colnames(data14)[20] <- "C_Pos"
colnames(data14)[6] <- "P-value"
write.table(data14, file = "/home/renshida/eqtm/data/data_cpg/32493484-cis.txt", sep = "\t", row.names = FALSE)




########33989148
# data15 <- read.table("33989148-cis.txt",header = TRUE,sep="")
# #data15 <- data15 %>% select(1)
# data15 <- merge(data15, 
#                 ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
#                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
# colnames(data15)[2] <- "CpG_Gene"
# colnames(data15)[3] <- "Feature"
# colnames(data15)[4] <- "Chr"
# colnames(data15)[5] <- "Pos"






###########30390659
data16 <- read.table("30390659-cis.txt",header = TRUE,sep = "")
#data16 <- data16 %>% select(1)
data16 <- merge(data16, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data16)[17] <- "CpG_Gene"
colnames(data16)[18] <- "Feature"
colnames(data16)[19] <- "C_Chr"
colnames(data16)[20] <- "C_Pos"
colnames(data16)[6] <- "P-value"
write.table(data16, file = "/home/renshida/eqtm/data/data_cpg/30390659-cis.txt", sep = "\t", row.names = FALSE)


################36803404
data17 <- read.table("36803404-top25.txt",header = TRUE,sep = "")
#data17 <- data17 %>% select(1)
data17 <- merge(data17, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data17)[17] <- "CpG_Gene"
colnames(data17)[18] <- "Feature"
colnames(data17)[19] <- "C_Chr"
colnames(data17)[20] <- "C_Pos"
colnames(data17)[6] <- "P-value"
str(data17$Beta)
data17$Beta <- gsub("−", "-", data17$Beta)
data17$Beta <- gsub("\\s+", "", data17$Beta)
data17$`P-value` <- as.numeric(data17$`P-value`)
data17$Beta <- as.numeric(data17$Beta)
data17$z <- abs(qnorm(data17$`P-value` / 2))
data17$SE <- abs(data17$Beta) / data17$z
data17 <- data17[,-21]
write.table(data17, file = "/home/renshida/eqtm/data/data_cpg/36803404-top25.txt", sep = "\t", row.names = FALSE)




##########32569636
data18 <- read.table("32569636-cis-top30-eva.txt",header = TRUE,sep="")
#data18 <- data18 %>% select(1)
data18 <- merge(data18, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data18)[17] <- "CpG_Gene"
colnames(data18)[18] <- "Feature"
colnames(data18)[19] <- "C_Chr"
colnames(data18)[20] <- "C_Pos"
colnames(data18)[6] <- "P-value"
str(data18$`P-value`)
data18$`P-value` <- gsub(" × 10", "e", data18$`P-value`)
data18$`P-value` <- as.numeric(data18$`P-value`)
data18$Beta <- as.numeric(data18$Beta)
data18$z <- abs(qnorm(data18$`P-value` / 2))
data18$SE <- abs(data18$Beta) / data18$z
data18 <- data18[,-21]
write.table(data18, file = "/home/renshida/eqtm/data/data_cpg/32569636-cis-top30-eva.txt", sep = "\t", row.names = FALSE)

data18.1 <- read.table("32569636-cis-top30-gse.txt",header = TRUE,sep="")
#data18.1 <- data18.1 %>% select(1)
data18.1 <- merge(data18.1, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data18.1)[17] <- "CpG_Gene"
colnames(data18.1)[18] <- "Feature"
colnames(data18.1)[19] <- "C_Chr"
colnames(data18.1)[20] <- "C_Pos"
colnames(data18.1)[6] <- "P-value"
str(data18.1$`P-value`)
data18.1$`P-value` <- gsub(" × 10", "e", data18.1$`P-value`)
data18.1$`P-value` <- as.numeric(data18.1$`P-value`)
data18.1$Beta <- as.numeric(data18.1$Beta)
data18.1$z <- abs(qnorm(data18.1$`P-value` / 2))
data18.1$SE <- abs(data18.1$Beta) / data18.1$z
data18.1 <- data18.1[,-21]
write.table(data18.1, file = "/home/renshida/eqtm/data/data_cpg/32569636-cis-top30-gse.txt", sep = "\t", row.names = FALSE)




#####
data14_parsed$h38_chr <- data14_parsed$Chromosome
data14_parsed$h38_pos <- data14_parsed$Position

# 使用 merge 函数进行匹配
merged_data <- NULL

# 保留需要的列
result <- NULL

######37563237
data19 <- read.table("37563237-cis-top10000.txt",header = TRUE,sep="")
#data18 <- data18 %>% select(1)
data19 <- merge(data19, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data19)[17] <- "CpG_Gene"
colnames(data19)[18] <- "Feature"
colnames(data19)[19] <- "C_Chr"
colnames(data19)[20] <- "C_Pos"
colnames(data19)[6] <- "P-value"
data19$`P-value` <- as.numeric(data19$`P-value`)
data19$Beta <- as.numeric(data19$Beta)
data19$z <- abs(qnorm(data19$`P-value` / 2))
data19$SE <- abs(data19$Beta) / data19$z
data19 <- data19[,-21]
write.table(data19, file = "/home/renshida/eqtm/data/data_cpg/37563237-cis-top10000.txt", sep = "\t", row.names = FALSE)



data19.1 <- read.table("37563237-trans-top10000.txt",header = TRUE,sep="")
#data18 <- data18 %>% select(1)
data19.1 <- merge(data19.1, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data19.1)[17] <- "CpG_Gene"
colnames(data19.1)[18] <- "Feature"
colnames(data19.1)[19] <- "C_Chr"
colnames(data19.1)[20] <- "C_Pos"
colnames(data19.1)[6] <- "P-value"
data19.1$`P-value` <- as.numeric(data19.1$`P-value`)
data19.1$Beta <- as.numeric(data19.1$Beta)
data19.1$z <- abs(qnorm(data19.1$`P-value` / 2))
data19.1$SE <- abs(data19.1$Beta) / data19.1$z
data19.1 <- data19.1[,-21]
write.table(data19.1, file = "/home/renshida/eqtm/data/data_cpg/37563237-trans-top10000.txt", sep = "\t", row.names = FALSE)


setwd("/home/renshida/eqtm/data/data_cpg")
#####32928872
data20 <- read.table("32928872-cis-cite.txt",header = TRUE,sep="")
colnames(data20)[6] <- "P-value"
write.table(data20, file = "/home/renshida/eqtm/data/data_cpg/32928872-cis-cite.txt", sep = "\t", row.names = FALSE)

data20.1 <- read.table("32928872-cis-region.txt",header = TRUE,sep="")
colnames(data20.1)[7] <- "P-value"
write.table(data20.1, file = "/home/renshida/eqtm/data/data_cpg/32928872-cis-region.txt", sep = "\t", row.names = FALSE)


#####31403346
data21 <- read.table("31403346-cis.txt",header = TRUE,sep="")
colnames(data21)[6] <- "P-value"
data21$`P-value` <- as.numeric(data21$`P-value`)
data21$Beta <- as.numeric(data21$Beta)
data21$z <- abs(qnorm(data21$`P-value` / 2))
data21$SE <- abs(data21$Beta) / data21$z
data21 <- data21[,-17]
write.table(data21, file = "/home/renshida/eqtm/data/data_cpg/31403346-cis.txt", sep = "\t", row.names = FALSE)


####37047575
data22 <- read.table("37047575-cis-DMR.txt",header = TRUE,sep="")
colnames(data22)[7] <- "P-value"
data22$`P-value` <- as.numeric(data22$`P-value`)
data22$Beta <- as.numeric(data22$Beta)
data22$z <- abs(qnorm(data22$`P-value` / 2))
data22$SE <- abs(data22$Beta) / data22$z
data22 <- data22[,-18]
write.table(data22, file = "/home/renshida/eqtm/data/data_cpg/37047575-cis-DMR.txt", sep = "\t", row.names = FALSE)



#####31578227
data23 <- read.table("31578227-cis-DMR.txt",header = TRUE,sep="")
colnames(data23)[7] <- "P-value"

data23$`P-value` <- as.numeric(data23$`P-value`)
data23$Beta <- as.numeric(data23$Beta)
data23$z <- abs(qnorm(data23$`P-value` / 2))
data23$SE <- abs(data23$Beta) / data23$z
data23 <- data23[,-18]
write.table(data23, file = "/home/renshida/eqtm/data/data_cpg/31578227-cis-DMR.txt", sep = "\t", row.names = FALSE)



######35232286
data24 <- read.table("35232286-DMR.txt",header = TRUE,sep="")
colnames(data24)[7] <- "P-value"
data24 <- data24[c(-1527,-1528),]
data24$`P-value` <- as.numeric(data24$`P-value`)
data24$Beta <- as.numeric(data24$Beta)
data24$z <- abs(qnorm(data24$`P-value` / 2))
data24$SE <- abs(data24$Beta) / data24$z
data24 <- data24[,-18]
write.table(data24, file = "/home/renshida/eqtm/data/data_cpg/35232286-DMR.txt", sep = "\t", row.names = FALSE)


















