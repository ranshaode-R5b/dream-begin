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
data1 <- data1 %>% select(1)
data1 <- merge(data1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data1)[2] <- "CpG_Gene"
colnames(data1)[3] <- "Feature"
colnames(data1)[4] <- "Chr"
colnames(data1)[5] <- "Pos"
head(data1)


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
data2 <- data2 %>% select(1)
data2 <- merge(data2, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2)[2] <- "CpG_Gene"
colnames(data2)[3] <- "Feature"
colnames(data2)[4] <- "Chr"
colnames(data2)[5] <- "Pos"


data2.1 <- read.table("37271320-trans1.txt", header = TRUE, sep = "")

# 仅保留第一列
data2.1 <- data2.1 %>% select(1)
data2.1 <- merge(data2.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2.1)[2] <- "CpG_Gene"
colnames(data2.1)[3] <- "Feature"
colnames(data2.1)[4] <- "Chr"
colnames(data2.1)[5] <- "Pos"


data2.2 <- read.table("37271320-trans2.txt", header = TRUE, sep = "")
# 仅保留第一列
data2.2 <- data2.2 %>% select(1)
data2.2 <- merge(data2.2, 
                 ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2.2)[2] <- "CpG_Gene"
colnames(data2.2)[3] <- "Feature"
colnames(data2.2)[4] <- "Chr"
colnames(data2.2)[5] <- "Pos"

data2.3 <- read.table("37271320-trans3.txt", header = TRUE, sep = "")
# 仅保留第一列
data2.3 <- data2.3 %>% select(1)
data2.3 <- merge(data2.3, 
                 ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data2.3)[2] <- "CpG_Gene"
colnames(data2.3)[3] <- "Feature"
colnames(data2.3)[4] <- "Chr"
colnames(data2.3)[5] <- "Pos"


########37233989
data3 <- read.table("37233989-cis-vc.txt", header = TRUE, sep = "")

# 仅保留第一列
data3 <- data3 %>% select(1)
data3 <- merge(data3, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data3)[2] <- "CpG_Gene"
colnames(data3)[3] <- "Feature"
colnames(data3)[4] <- "Chr"
colnames(data3)[5] <- "Pos"

data3.1 <- read.table("37233989-cis-ve.txt", header = TRUE, sep = "")
# 仅保留第一列
data3.1 <- data3.1 %>% select(1)
data3.1 <- merge(data3.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data3.1)[2] <- "CpG_Gene"
colnames(data3.1)[3] <- "Feature"
colnames(data3.1)[4] <- "Chr"
colnames(data3.1)[5] <- "Pos"




#########39393618
data4 <- read.table("39393618-cis.txt",header = TRUE,sep="")
data4 <- data4 %>% select(1)
data4 <- merge(data4, 
                 ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data4)[2] <- "CpG_Gene"
colnames(data4)[3] <- "Feature"
colnames(data4)[4] <- "Chr"
colnames(data4)[5] <- "Pos"




#######31076557
data5 <- read.table("31076557-cis.txt",header = TRUE,sep="")
data5 <- data5 %>% select(1)
data5 <- merge(data5, 
               ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data5)[2] <- "CpG_Gene"
colnames(data5)[3] <- "Feature"
colnames(data5)[4] <- "Chr"
colnames(data5)[5] <- "Pos"



########35710981
data6 <- read.table("35710981-cis.txt",header = TRUE,sep="")
data6 <- data6 %>% select(1)
data6 <- merge(data6, 
               ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data6)[2] <- "CpG_Gene"
colnames(data6)[3] <- "Feature"
colnames(data6)[4] <- "Chr"
colnames(data6)[5] <- "Pos"


#######31791327
data7 <- read.table("31791327-cis-LL.txt",header = TRUE,sep="")

data7 <- data7 %>% select(1)
data7 <- merge(data7, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7)[2] <- "CpG_Gene"
colnames(data7)[3] <- "Feature"
colnames(data7)[4] <- "Chr"
colnames(data7)[5] <- "Pos"


data7.1 <- read.table("31791327-cis-LLS.txt",header = TRUE,sep="")

data7.1 <- data7.1 %>% select(1)
data7.1 <- merge(data7.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7.1)[2] <- "CpG_Gene"
colnames(data7.1)[3] <- "Feature"
colnames(data7.1)[4] <- "Chr"
colnames(data7.1)[5] <- "Pos"


data7.2 <- read.table("31791327-cis-NTR.txt",header = TRUE,sep="")

data7.2 <- data7.2 %>% select(1)
data7.2 <- merge(data7.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7.2)[2] <- "CpG_Gene"
colnames(data7.2)[3] <- "Feature"
colnames(data7.2)[4] <- "Chr"
colnames(data7.2)[5] <- "Pos"



data7.3 <- read.table("31791327-cis-RS.txt",header = TRUE,sep="")

data7.3 <- data7.3 %>% select(1)
data7.3 <- merge(data7.3, 
                 ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                 by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data7.3)[2] <- "CpG_Gene"
colnames(data7.3)[3] <- "Feature"
colnames(data7.3)[4] <- "Chr"
colnames(data7.3)[5] <- "Pos"


###############35302492
data8 <- read.table("35302492-cis-adj.txt",header = TRUE,sep="")

data8 <- data8 %>% select(1)
data8 <- merge(data8, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data8)[2] <- "CpG_Gene"
colnames(data8)[3] <- "Feature"
colnames(data8)[4] <- "Chr"
colnames(data8)[5] <- "Pos"


data8.1 <- read.table("35302492-cis-unadj.txt",header = TRUE,sep="")

data8.1 <- data8.1 %>% select(1)
data8.1 <- merge(data8.1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data8.1)[2] <- "CpG_Gene"
colnames(data8.1)[3] <- "Feature"
colnames(data8.1)[4] <- "Chr"
colnames(data8.1)[5] <- "Pos"




######37598461
data9 <- read.table("37598461-cis.txt",header = TRUE,sep="")
data9 <- data9 %>% select(1)
data9 <- merge(data9, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data9)[2] <- "CpG_Gene"
colnames(data9)[3] <- "Feature"
colnames(data9)[4] <- "Chr"
colnames(data9)[5] <- "Pos"




#######33092652
data10 <- read.table("33092652-cis-trans-kora.txt",header = TRUE,sep="")
colnames(data10)[1] <- "CpG"
data10 <- data10 %>% select(1)
data10 <- merge(data10, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data10)[2] <- "CpG_Gene"
colnames(data10)[3] <- "Feature"
colnames(data10)[4] <- "Chr"
colnames(data10)[5] <- "Pos"


data10.1 <- read.table("33092652-cis-trans-rot.txt",header = TRUE,sep="")
colnames(data10.1)[1] <- "CpG"
data10.1 <- data10.1 %>% select(1)
data10.1 <- merge(data10.1, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data10.1)[2] <- "CpG_Gene"
colnames(data10.1)[3] <- "Feature"
colnames(data10.1)[4] <- "Chr"
colnames(data10.1)[5] <- "Pos"




#######29914364
data11 <- read.table("29914364-GTP.txt",header = TRUE,sep="")
data11 <- data11 %>% select(1)
data11 <- merge(data11, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data11)[2] <- "CpG_Gene"
colnames(data11)[3] <- "Feature"
colnames(data11)[4] <- "Chr"
colnames(data11)[5] <- "Pos"


data11.1 <- read.table("29914364-MESA.txt",header = TRUE,sep="")
data11.1 <- data11.1 %>% select(1)
data11.1 <- merge(data11.1, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data11.1)[2] <- "CpG_Gene"
colnames(data11.1)[3] <- "Feature"
colnames(data11.1)[4] <- "Chr"
colnames(data11.1)[5] <- "Pos"




#############35426765
data12 <- read.table("35426765-cis.txt",header = TRUE,sep = "")
colnames(data12)[1] <- "CpG"
data12 <- data12 %>% select(1)
data12 <- merge(data12, 
               ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data12)[2] <- "CpG_Gene"
colnames(data12)[3] <- "Feature"
colnames(data12)[4] <- "Chr"
colnames(data12)[5] <- "Pos"



###########38438351
data13 <- read.table("38438351-cis.txt",header = TRUE,sep = "")
colnames(data13)[1] <- "CpG"
data13 <- data13 %>% select(1)
data13 <- merge(data13, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data13)[2] <- "CpG_Gene"
colnames(data13)[3] <- "Feature"
colnames(data13)[4] <- "Chr"
colnames(data13)[5] <- "Pos"




######32493484
data14 <- read.table("32493484-cis.txt",header = TRUE,sep="")
colnames(data14)[1] <- "CpG"
data14 <- data14 %>% select(1)
data14 <- merge(data14, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data14)[2] <- "CpG_Gene"
colnames(data14)[3] <- "Feature"
colnames(data14)[4] <- "Chr"
colnames(data14)[5] <- "Pos"





########33989148
data15 <- read.table("33989148-cis.txt",header = TRUE,sep="")
data15 <- data15 %>% select(1)
data15 <- merge(data15, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data15)[2] <- "CpG_Gene"
colnames(data15)[3] <- "Feature"
colnames(data15)[4] <- "Chr"
colnames(data15)[5] <- "Pos"






###########30390659
data16 <- read.table("30390659-cis.txt",header = TRUE,sep = "")
data16 <- data16 %>% select(1)
data16 <- merge(data16, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data16)[2] <- "CpG_Gene"
colnames(data16)[3] <- "Feature"
colnames(data16)[4] <- "Chr"
colnames(data16)[5] <- "Pos"



################36803404
data17 <- read.table("36803404-top25.txt",header = TRUE,sep = "")
data17 <- data17 %>% select(1)
data17 <- merge(data17, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data17)[2] <- "CpG_Gene"
colnames(data17)[3] <- "Feature"
colnames(data17)[4] <- "Chr"
colnames(data17)[5] <- "Pos"





##########32569636
data18 <- read.table("32569636-cis-top30-eva.txt",header = TRUE,sep="")
data18 <- data18 %>% select(1)
data18 <- merge(data18, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data18)[2] <- "CpG_Gene"
colnames(data18)[3] <- "Feature"
colnames(data18)[4] <- "Chr"
colnames(data18)[5] <- "Pos"


data18.1 <- read.table("32569636-cis-top30-gse.txt",header = TRUE,sep="")
data18.1 <- data18.1 %>% select(1)
data18.1 <- merge(data18.1, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "Relation_to_UCSC_CpG_Island", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data18.1)[2] <- "CpG_Gene"
colnames(data18.1)[3] <- "Feature"
colnames(data18.1)[4] <- "Chr"
colnames(data18.1)[5] <- "Pos"





#####
data14_parsed$h38_chr <- data14_parsed$Chromosome
data14_parsed$h38_pos <- data14_parsed$Position

# 使用 merge 函数进行匹配
merged_data <- NULL

# 保留需要的列
result <- NULL

























