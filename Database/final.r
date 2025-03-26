setwd("/home/renshida/eqtm/data/data_cpg/")
######37271320
data1 <- read.table("37271320-all.txt",header = T,sep = "")
data1 <- data1 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
colnames(data1)
colnames(ann_450k)
a <- c(2,1,17,19,20,21,18,3:16)
data1 <- data1[,a]
data1 <- data1 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data1 <- data1[,-c(5,6)]
colnames(data1)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data1, file = "/data1/renshida/00_final_eqtm/37271320-all.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)
range(data1$`P-value`)
range(data1$FDR)



#######29914364
data2 <- read.table("29914364-GTP.txt",header = T,sep = "")
data2$Date <- "2018.06"
range(data2$P.value)
data2 <- data2 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data2 <- data2[,a]
data2 <- data2 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data2 <- data2[,-c(5,6)]
colnames(data2)[6] <- "Ref_Gene_Group"
colnames(data2)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data2, file = "/data1/renshida/00_final_eqtm/29914364-GTP.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)
colnames(data2)

data2.1 <- read.table("29914364-MESA.txt",header = T,sep = "")
data2.1$Date <- "2018.06"
range(data2.1$P.value)
data2.1 <- data2.1 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data2.1 <- data2.1[,a]
data2.1 <- data2.1 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data2.1 <- data2.1[,-c(5,6)]
colnames(data2.1)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data2.1, file = "/data1/renshida/00_final_eqtm/29914364-MESA.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


#####35426765
data3 <- read.table("35426765-cis.txt",header = T,sep = "")
range(data3$P.value)
range(data3$FDR)
data3 <- data3 %>%
  left_join(ann_EPIC %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data3 <- data3[,a]
data3 <- data3 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data3 <- data3[,-c(5,6)]
colnames(data3)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data3, file = "/data1/renshida/00_final_eqtm/35426765-cis.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


####35302492
data4 <- read_tsv("/home/renshida/eqtm/rawdata/human/35302492-cis-unadj.txt")
unique_gene_count <- data4 %>%
  filter(TC_gene == "" | is.na(TC_gene)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(TC) %>%                                 # 去重 Gene 列
  nrow() 
unique_gene_count
length(unique(data4$TC))
data4$TC_gene <- sub(";.*", "", data4$TC_gene)
#HTA2.0 <- read.table("/home/renshida/eqtm/rawdata/human/HTA2.0.txt",header = TRUE, sep = "\t")
#head(HTA2.0)
data4 <- data4 %>%
  left_join(HTA2.0 %>% select(ID, X), by = c("TC" = "ID"))
data4$X[data4$X == "---"] <- NA
unique_gene_count <- data4 %>%
  filter(X == "" | is.na(X)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(TC) %>%                                 # 去重 Gene 列
  nrow()
unique_gene_count
length(unique(data4$TC))
data4 <- data4 %>%
  mutate(CombinedGeneName = coalesce(TC_gene, X))
unique_gene_count <- data4 %>%
  filter(CombinedGeneName == "" | is.na(CombinedGeneName)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(TC) %>%                                 # 去重 Gene 列
  nrow()
unique_gene_count
length(unique(data4$TC))
colnames(data4)[1] <- "CpG"
colnames(data4)[3] <- "Beta"
colnames(data4)[4] <- "SE"
colnames(data4)[5] <- "P-value"
data4 <- data4[,c(1,18,3,4,5)]
colnames(data4)[2] <- "Gene"
data4$Species <- "Human"
data4$Tissue <- "Whole blood"
data4$Disease <- "Autosomal"
data4$Article <- "Identification of autosomal cis expression quantitative trait methylation (cis eQTMs) in children's blood"
data4$Journal <- "Elife"
data4$Date <- "2022.03"
data4$PMID <- "35302492"
data4$Cis_Trans_eqtm <- "Cis"
data4$Method <- "Linear regression model"
data4$Sample <- "832"
data4$FDR <- NA
a <- c(1:2,13,3:12,14:16)
data4 <- data4[,a]
a <- c(1:6,16,7:15)
data4 <- data4[,a]
data4 <- merge(data4, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
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
data4$`P-value` <- as.numeric(data4$`P-value`)
data4 <- data4 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data4 <- data4[,a]
data4 <- data4 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data4 <- data4[,-c(5,6)]
colnames(data4)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
data4 <- data4[!is.na(data4$`Gene Symbol`), ]
write.table(data4, file = "/data1/renshida/00_final_eqtm/35302492-cis-unadj.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


#######35710981
data5 <- read.table("35710981-cis.txt",header = T,sep = "")
range(data5$P.value)
range(data5$FDR)
data5 <- data5 %>%
  left_join(ann_EPIC %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data5 <- data5[,a]
data5 <- data5 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data5 <- data5[,-c(5,6)]
colnames(data5)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data5, file = "/data1/renshida/00_final_eqtm/35710981-cis.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


######31076557
data6 <- read_tsv("/home/renshida/eqtm/rawdata/human/31076557-cis.txt")
data6 <- data6[!is.na(data6$pv), ]
range(data6$pv)
range(data6$qv_bh)
data6 <- data6[,c(6,3,7,9,10)]
data6 <- data6[!is.na(data6$pv), ]
colnames(data6)[1] <- "CpG"
colnames(data6)[2] <- "Gene"
colnames(data6)[3] <- "Beta"
colnames(data6)[4] <- "P-value"
colnames(data6)[5] <- "FDR"
data6 <- data6 %>% filter(FDR < 0.01)
data6$Species <- "Human"
data6$Tissue <- "Skeletal muscle"
data6$Disease <- "Normal"
data6$Article <- "Integrative analysis of gene expression, DNA methylation, physiological traits, and genetic variation in human skeletal muscle"
data6$Journal <- "Proceedings of the National Academy of Sciences of the United States of America"
data6$Date <- "2019.05"
data6$PMID <- "31076557"
data6$Cis_Trans_eqtm <- "Cis"
data6$Method <- "Linear regression model"
data6$Sample <- 318
data6$SE <- NA
a <- c(1,2,13,3:12,14,15,16)
data6 <- data6[,a]
a <- c(1:4,16,5:15)
data6 <- data6[,a]
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
range(data6$`P-value`)
range(data6$FDR)
data6 <- data6 %>%
  left_join(ann_EPIC %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data6 <- data6[,a]
data6 <- data6 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data6 <- data6[,-c(5,6)]
colnames(data6)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
gtf <- read.table("/home/renshida/eqtm/rawdata/human/extracted_columns_hg19.txt")
colnames(gtf)[1] <- "ens"
colnames(gtf)[2] <- "sym"
head(gtf)
data6$`Gene Symbol` <- sub("\\..*", "", data6$`Gene Symbol`)
gtf$ens <- sub("\\..*", "", gtf$ens)
data6 <- data6 %>%
  left_join(gtf, by = c(`Gene Symbol` = "ens"))
a <- c(21,2:20)
data6 <- data6[,a]
colnames(data6)[1] <- "Gene Symbol"
#data6 <- data6[!is.na(data6$Gene), ]
has_na <- any(is.na(data6$`Gene Symbol`))
na_count <- sum(is.na(data6$`Gene Symbol`))
print(has_na)
print(na_count)
write.table(data6, file = "/data1/renshida/00_final_eqtm/31076557-cis.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


####39393618
data7 <- read.table("39393618-cis.txt",header = T,sep = "")
range(data7$P.value)
range(data7$FDR)
data7 <- data7 %>%
  left_join(ann_EPIC %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data7 <- data7[,a]
data7 <- data7 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data7 <- data7[,-c(5,6)]
colnames(data7)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data7, file = "/data1/renshida/00_final_eqtm/39393618-cis.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)



#######33092652
file_path <- "/home/renshida/eqtm/rawdata/human/33092652-cis-trans.xlsx"
sheets <- excel_sheets(file_path)
data_list <- lapply(sheets, function(sheet) {
  read_excel(file_path, sheet = sheet)
})
# 初始化一个空列表来存储处理过的数据框
# 遍历每个工作表
processed_data_list <- list()

# 遍历每个工作表
for (i in seq_along(data_list)) {
  data <- data_list[[i]]
  
  # 获取第7列的列名
  gene_column_name <- colnames(data)[7]
  
  # 删除第一行和第二行
  data <- data[-c(1, 2), ]
  
  # 将第7列的列名作为新列添加到每一行
  data$Gene <- gene_column_name
  
  # 将处理过的数据框添加到列表中
  processed_data_list[[i]] <- data
}
data8 <- bind_rows(processed_data_list)
colnames(data8)[1] <- "CpG"
colnames(data8)[2] <- "Beta"
colnames(data8)[3] <- "P-value"
a <- c(1,8,2:7,9:ncol(data8))
data8 <- data8[, a]
data8 <- data8[,c(1,2,3,4)]
data8$SE <- NA
data8$FDR <- NA
data8$Species <- "Human"
data8$Tissue <- "Whole blood"
data8$Disease <- "Differences between smokers and never smokers"
data8$Article <- "Smoking-related changes in DNA methylation and gene expression are associated with cardio-metabolic traits"
data8$Journal <- "Clinical Epigenetics"
data8$Date <- "2020.10"
data8$PMID <- "33092652"
data8$Cis_Trans_eqtm <- NA
data8$Method <- "Linear mixed model "
data8$Sample <- "716"
#na_columns <- apply(combined_data, 2, function(x) any(is.na(x)))
#print(na_columns)
a <- c(1:3,5,4,6:ncol(data8))
data8 <- data8[, a]
data8$'P-value' <- as.numeric(data8$`P-value`)
pvalue_range <- range(data8$'P-value', na.rm = TRUE)
print(pvalue_range)
a <- c(1:2,14,3:13,15,16)
data8 <- data8[,a]
data8 <- merge(data8, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data8)[17] <- "CpG_Gene"
colnames(data8)[18] <- "Feature"
colnames(data8)[19] <- "C_Chr"
colnames(data8)[20] <- "C_Pos"
colnames(data8)[6] <- "P-value"
data8$`P-value` <- as.numeric(data8$`P-value`)
data8$Beta <- as.numeric(data8$Beta)
data8$z <- abs(qnorm(data8$`P-value` / 2))
data8$SE <- abs(data8$Beta) / data8$z
data8 <- data8[,-21]


probe_ids <- data8$Gene
# 利用 mapIds 函数将 Illumina 探针 ID 转换为 Gene Symbol
gene_symbols <- mapIds(illuminaHumanv4.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# 将转换结果添加为新的列
data8$GeneSymbol <- gene_symbols
has_na <- any(is.na(data8$GeneSymbol))
na_count <- sum(is.na(data8$GeneSymbol))
print(has_na)
print(na_count)
unique_gene_count <- data8 %>%
  filter(GeneSymbol == "" | is.na(GeneSymbol)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(Gene) %>%                                 # 去重 Gene 列
  nrow() 
unique_gene_count
length(unique(data8$Gene))
a <- c(1,21,3:20)
data8 <- data8[,a]
colnames(data8)[2] <- "Gene"
length(unique(data8$Gene))
##data8 <- data8[!is.na(data8$Gene), ]


data8 <- data8 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data8 <- data8[,a]
data8 <- data8 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data8 <- data8[,-c(5,6)]
colnames(data8)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data8, file = "/data1/renshida/00_final_eqtm/33092652-rot.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)
has_na <- any(is.na(data8$`Gene Symbol`))
na_count <- sum(is.na(data8$`Gene Symbol`))
print(has_na)
print(na_count)




combined_data1 <- bind_rows(processed_data_list)
colnames(combined_data1)[1] <- "CpG"
colnames(combined_data1)[4] <- "Beta"
colnames(combined_data1)[5] <- "P-value"
a <- c(1,8,2:7,9:ncol(combined_data1))
combined_data1 <- combined_data1[, a]
combined_data1 <- combined_data1[,c(1,2,5,6)]
#na_columns <- apply(combined_data1, 2, function(x) any(is.na(x)))
#print(na_columns)
combined_data1 <- na.omit(combined_data1)
combined_data1$SE <- NA
combined_data1$FDR <- NA
combined_data1$Species <- "Human"
combined_data1$Tissue <- "Whole blood"
combined_data1$Disease <- "Differences between smokers and never smokers"
combined_data1$Article <- "Smoking-related changes in DNA methylation and gene expression are associated with cardio-metabolic traits"
combined_data1$Journal <- "Clinical Epigenetics"
combined_data1$Date <- "2020.10"
combined_data1$PMID <- "33092652"
combined_data1$Cis_Trans_eqtm <- NA
combined_data1$Method <- "Linear mixed model "
combined_data1$Sample <- "687"
#na_columns <- apply(combined_data, 2, function(x) any(is.na(x)))
#print(na_columns)
a <- c(1:3,5,4,6:ncol(combined_data1))
combined_data1 <- combined_data1[, a]
combined_data1$'P-value' <- as.numeric(combined_data1$`P-value`)
pvalue_range <- range(combined_data1$'P-value', na.rm = TRUE)
print(pvalue_range)
a <- c(1:2,14,3:13,15,16)
combined_data1 <- combined_data1[,a]

combined_data1 <- merge(combined_data1, 
               ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
               by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(combined_data1)[17] <- "CpG_Gene"
colnames(combined_data1)[18] <- "Feature"
colnames(combined_data1)[19] <- "C_Chr"
colnames(combined_data1)[20] <- "C_Pos"
colnames(combined_data1)[6] <- "P-value"
combined_data1$`P-value` <- as.numeric(combined_data1$`P-value`)
combined_data1$Beta <- as.numeric(combined_data1$Beta)
combined_data1$z <- abs(qnorm(combined_data1$`P-value` / 2))
combined_data1$SE <- abs(combined_data1$Beta) / combined_data1$z
combined_data1 <- combined_data1[,-21]


probe_ids <- combined_data1$Gene
# 利用 mapIds 函数将 Illumina 探针 ID 转换为 Gene Symbol
gene_symbols <- mapIds(illuminaHumanv4.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# 将转换结果添加为新的列
combined_data1$GeneSymbol <- gene_symbols
has_na <- any(is.na(combined_data1$GeneSymbol))
na_count <- sum(is.na(combined_data1$GeneSymbol))
print(has_na)
print(na_count)

unique_gene_count <- combined_data1 %>%
  filter(GeneSymbol == "" | is.na(GeneSymbol)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(Gene) %>%                                 # 去重 Gene 列
  nrow() 
unique_gene_count
length(unique(combined_data1$Gene))
a <- c(1,21,3:20)
combined_data1 <- combined_data1[,a]
colnames(combined_data1)[2] <- "Gene"
##data8 <- data8[!is.na(data8$Gene), ]


combined_data1 <- combined_data1 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
combined_data1 <- combined_data1[,a]
combined_data1 <- combined_data1 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
combined_data1 <- combined_data1[,-c(5,6)]
colnames(combined_data1)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(combined_data1, file = "/data1/renshida/00_final_eqtm/33092652-kora.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)
has_na <- any(is.na(combined_data1$`Gene Symbol`))
na_count <- sum(is.na(combined_data1$`Gene Symbol`))
print(has_na)
print(na_count)




########37563237
data9.1 <- read.table("37563237-cis-top10000.txt",header = T,sep = "")
data9.2 <- read.table("37563237-trans-top10000.txt",header = T,sep = "")
data9 <- rbind(data9.1,data9.2)
data9$P.value <- as.numeric(data9$P.value)
range(data9.2$P.value)
range(data9$FDR)
data9 <- data9 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data9 <- data9[,a]
data9 <- data9 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data9 <- data9[,-c(5,6)]
colnames(data9)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data9, file = "/data1/renshida/00_final_eqtm/37563237-ALL.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)





#######38438351
data10 <- read.table("/home/renshida/eqtm/rawdata/human/38438351-cis.txt",header = T,sep = "")
data10 <- data10[,c(2,1,4,5,3,9)]
colnames(data10)[1] <- "CpG"
colnames(data10)[2] <- "Gene"
colnames(data10)[3] <- "Beta"
colnames(data10)[4] <- "SE"
colnames(data10)[5] <- "P-value"
colnames(data10)[6] <- "FDR"
data10$Species <- "Human"
data10$Tissue <- "Retina"
data10$Disease <- "Normal"
data10$Article <- "QTL mapping of human retina DNA methylation identifies 87 gene-epigenome interactions in age-related macular degeneration"
data10$Journal <- "Nature Communications"
data10$Date <- "2024.03"
data10$PMID <- "38438351"
data10$Cis_Trans_eqtm <- "Cis"
data10$Method <- "Linear regression model"
data10$Sample <- "152"
a <- c(1:2,14,3:13,15,16)
data10 <- data10[,a]
data10$Gene <- sub("\\..*", "", data10$Gene)
#gtf <- read.table("/home/renshida/eqtm/rawdata/human/gene_id_name.txt")
colnames(gtf)[1] <- "ens"
colnames(gtf)[2] <- "sym"
gtf$ens <- sub("\\..*", "", gtf$ens)
data10 <- data10 %>%
  left_join(gtf, by = c("Gene" = "ens"))
a <- c(1,17,3:16)
data10 <- data10[,a]
colnames(data10)[2] <- "GeneSymbol"
has_na <- any(is.na(data10$`GeneSymbol`))
na_count <- sum(is.na(data10$`GeneSymbol`))
print(has_na)
print(na_count)
# unique_gene_count <- data10 %>%
#   filter(GeneSymbol == "" | is.na(GeneSymbol)) %>%  # 过滤 Genesymbol 为空或 NA 的行
#   distinct(Gene) %>%                                 # 去重 Gene 列
#   nrow() 
# unique_gene_count
# length(unique(data10$Gene))
data10 <- data10[!is.na(data10$`GeneSymbol`), ]
data10$`P-value` <- as.numeric(data10$`P-value`)
data10$FDR <- as.numeric(data10$FDR)
range(data10$`P-value`)
range(data10$FDR)
data10 <- merge(data10, 
                ann_EPIC[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "CHR_hg38", "Start_hg38")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data10)[17] <- "CpG_Gene"
colnames(data10)[18] <- "Feature"
colnames(data10)[19] <- "C_Chr"
colnames(data10)[20] <- "C_Pos"
colnames(data10)[6] <- "P-value"
data10 <- data10 %>%
  left_join(ann_EPIC %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data10 <- data10[,a]
data10 <- data10 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(GeneSymbol, CpG, CpG_Gene, CpG_position, everything())
data10 <- data10[,-c(5,6)]
colnames(data10)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data10, file = "/data1/renshida/00_final_eqtm/38438351-cis.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


#####35232286
data11 <- read.table("35232286-DMR.txt",header = T,sep = "")
data11$P.value <- as.numeric(data11$P.value)
data11$Date <- "2022.04"
a <- c(3,1,2,4:17)
data11 <- data11[,a]
colnames(data11)[1] <- "Gene Symbol"
colnames(data11)[4] <- "Cis/Trans"
colnames(data11)[7] <- "P-value"
write.table(data11, file = "/data1/renshida/00_final_eqtm/35232286-DMR.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)
colnames(data11)



#####37233989
data12 <- read.table("37233989-cis-vc.txt",header = T,sep = "")
data12$P.value <- as.numeric(data12$P.value)
range(data12$P.value)
range(data12$FDR)
data12 <- data12 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data12 <- data12[,a]
data12 <- data12 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data12 <- data12[,-c(5,6)]
colnames(data12)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data12, file = "/data1/renshida/00_final_eqtm/37233989-cis-vc.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)



####31578227
data13 <- read.table("31578227-cis-DMR.txt",header = T,sep = "")
data13$P.value <- as.numeric(data13$P.value)
a <- c(3,1,2,4:17)
data13 <- data13[,a]
colnames(data13)[1] <- "Gene Symbol"
colnames(data13)[4] <- "Cis/Trans"
colnames(data13)[7] <- "P-value"
write.table(data13, file = "/data1/renshida/00_final_eqtm/31578227-cis-DMR.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)




#####31282290
data14 <- read.table("31282290-cis-mir.txt",header = T,sep = "")
data14$P.value <- as.numeric(data14$P.value)
range(data14$P.value)
range(data14$FDR)
data14 <- data14 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data14 <- data14[,a]
data14 <- data14 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data14 <- data14[,-c(5,6)]
colnames(data14)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data14, file = "/data1/renshida/00_final_eqtm/31282290-cis-mir.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)



#####30390659
data15 <- read_excel("/home/renshida/eqtm/rawdata/human/30390659-cis.xlsx")
colnames(data15) <- data15[1,]
data15 <- data15[-1,]
unique_gene_count <- data15 %>%
  filter(Genename == "" | is.na(Genename)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(Probe) %>%                                 # 去重 Gene 列
  nrow() 
unique_gene_count
length(unique(data15$Probe))
HU133 <- read.table("/home/renshida/eqtm/rawdata/human/HU133.txt",sep="\t")
data15 <- data15 %>%
  left_join(HU133 %>% select(V1,V2), by = c("Probe" = "V1"))
data15 <- data15 %>%
  mutate(CombinedGeneName = coalesce(Genename, V2))
unique_gene_count <- data15 %>%
  filter(CombinedGeneName == "" | is.na(CombinedGeneName)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(Probe) %>%                                 # 去重 Gene 列
  nrow()
unique_gene_count
length(unique(data15$Probe))
data15 <- data15[,c(1,13,2,3,4,5)]
colnames(data15)[1] <- "CpG"
colnames(data15)[2] <- "Gene"
colnames(data15)[3] <- "Beta"
colnames(data15)[4] <- "SE"
colnames(data15)[5] <- "P-value"
colnames(data15)[6] <- "FDR"
data15$Species <- "Human"
data15$Tissue <- "Lung"
data15$Disease <- "Differences between smokers and never smokers"
data15$Article <- "From blood to lung tissue: effect of cigarette smoke on DNA methylation and lung function"
data15$Journal <- "Respiratory Research"
data15$Date <- "2018.11"
data15$PMID <- "30390659"
data15$Cis_Trans_eqtm <- "Cis"
data15$Method <- "Linear regression model"
data15$Sample <- 36
a <- c(1:2,14,3:13,15:16)
data15 <- data15[,a]
data15 <- data15 %>%
  mutate(Gene = na_if(Gene, "")) %>%
  filter(!is.na(Gene))
data15 <- merge(data15, 
                        ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                        by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data15)[17] <- "CpG_Gene"
colnames(data15)[18] <- "Feature"
colnames(data15)[19] <- "C_Chr"
colnames(data15)[20] <- "C_Pos"
colnames(data15)[6] <- "P-value"
data15 <- data15 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data15 <- data15[,a]
data15 <- data15 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data15 <- data15[,-c(5,6)]
colnames(data15)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data15, file = "/data1/renshida/00_final_eqtm/30390659-cis.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


####32928872
data16 <- read.table("32928872-cis-cite.txt",header = T,sep="")
data16$P.value <- as.numeric(data16$P.value)
a <- c(2,1,3:16)
data16 <- data16[,a]
colnames(data16)[1] <- "Gene Symbol"
colnames(data16)[2] <- "CpG site"
colnames(data16)[3] <- "Cis/Trans"
colnames(data16)[6] <- "P-value"
write.table(data16, file = "/data1/renshida/00_final_eqtm/32928872-cis-cite.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)



#######################合并不同版本的GTF
gtf1 <- read.table("/home/renshida/eqtm/rawdata/human/hg19.ensembol.gtf")
gtf2 <- read.table("/home/renshida/eqtm/rawdata/human/hg38.ensemble.gtf")
gtf3 <- read.table("/home/renshida/eqtm/rawdata/human/hg38.genecode.gtf")
gtf1$V1 <- sub("\\..*", "", gtf1$V1)
gtf2$V1 <- sub("\\..*", "", gtf2$V1)
gtf3$V1 <- sub("\\..*", "", gtf3$V1)
head(gtf1)
head(gtf2)

combined_gtf <- bind_rows(gtf1, gtf2)
combined_gtf <- combined_gtf %>%
  distinct(V1, .keep_all = TRUE)

combined_gtf <- bind_rows(combined_gtf,gtf3)
combined_gtf <- combined_gtf %>%
  distinct(V1, .keep_all = TRUE)
gtf <- combined_gtf


######31791327
data17 <- read_excel("/home/renshida/eqtm/rawdata/human/31791327-cis.xlsx")
data17 <- data17[-c(806:811),]
colnames(data17) <- data17[2,]
data17 <- data17[-c(1,2),]
data17.1 <- data17[,c(1,3,7,8,9)]
colnames(data17.1)[1] <- "CpG"
colnames(data17.1)[2] <- "Gene"
colnames(data17.1)[3] <- "Beta"
colnames(data17.1)[4] <- "SE"
colnames(data17.1)[5] <- "P-value"
data17.1$FDR <- NA
data17.1$Species <- "Human"
data17.1$Tissue <- "Whole blood"
data17.1$Disease <- "Never smokers"
data17.1$Article <- "DNA methylation is associated with lung function in never smokers"
data17.1$Journal <- "Respiratory Research"
data17.1$Date <- "2019.12"
data17.1$PMID <- "31791327"
data17.1$Cis_Trans_eqtm <- "Cis"
data17.1$Method <- "Robust linear regression"
data17.1$Sample <- "166"
a <- c(1,2,14,3:13,15,16)
data17.1 <- data17.1[,a]
data17.1 <- data17.1 %>%
  left_join(gtf, by = c("Gene" = "ens"))
a <- c(1,17,3:16)
data17.1 <- data17.1[,a]
colnames(data17.1)[2] <- "GeneSymbol"
has_na <- any(is.na(data17.1$`GeneSymbol`))
na_count <- sum(is.na(data17.1$`GeneSymbol`))
print(has_na)
print(na_count)
data17.1 <- merge(data17.1, 
                  ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                  by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data17.1)[17] <- "CpG_Gene"
colnames(data17.1)[18] <- "Feature"
colnames(data17.1)[19] <- "C_Chr"
colnames(data17.1)[20] <- "C_Pos"
colnames(data17.1)[6] <- "P-value"
data17.1 <- data17.1 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data17.1 <- data17.1[,a]
data17.1 <- data17.1 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(GeneSymbol, CpG, CpG_Gene, CpG_position, everything())
data17.1 <- data17.1[,-c(5,6)]
colnames(data17.1)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data17.1, file = "/data1/renshida/00_final_eqtm//31791327-cis-LL.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)



data17.2 <- data17[,c(1,3,10,11,12)]
colnames(data17.2)[1] <- "CpG"
colnames(data17.2)[2] <- "Gene"
colnames(data17.2)[3] <- "Beta"
colnames(data17.2)[4] <- "SE"
colnames(data17.2)[5] <- "P-value"
data17.2$FDR <- NA
data17.2$Species <- "Human"
data17.2$Tissue <- "Whole blood"
data17.2$Disease <- "Never smokers"
data17.2$Article <- "DNA methylation is associated with lung function in never smokers"
data17.2$Journal <- "Respiratory Research"
data17.2$Date <- "2019.12"
data17.2$PMID <- "31791327"
data17.2$Cis_Trans_eqtm <- "Cis"
data17.2$Method <- "Robust linear regression"
data17.2$Sample <- "906"
a <- c(1,2,14,3:13,15,16)
data17.2 <- data17.2[,a]
data17.2 <- data17.2 %>%
  left_join(gtf, by = c("Gene" = "ens"))
a <- c(1,17,3:16)
data17.2 <- data17.2[,a]
colnames(data17.2)[2] <- "GeneSymbol"
has_na <- any(is.na(data17.2$`GeneSymbol`))
na_count <- sum(is.na(data17.2$`GeneSymbol`))
print(has_na)
print(na_count)
data17.2 <- merge(data17.2, 
                  ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                  by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data17.2)[17] <- "CpG_Gene"
colnames(data17.2)[18] <- "Feature"
colnames(data17.2)[19] <- "C_Chr"
colnames(data17.2)[20] <- "C_Pos"
colnames(data17.2)[6] <- "P-value"
data17.2 <- data17.2 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data17.2 <- data17.2[,a]
data17.2 <- data17.2 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(GeneSymbol, CpG, CpG_Gene, CpG_position, everything())
data17.2 <- data17.2[,-c(5,6)]
colnames(data17.2)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data17.2, file = "/data1/renshida/00_final_eqtm//31791327-cis-LLS.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)



data17.3 <- data17[,c(1,3,13,14,15)]
colnames(data17.3)[1] <- "CpG"
colnames(data17.3)[2] <- "Gene"
colnames(data17.3)[3] <- "Beta"
colnames(data17.3)[4] <- "SE"
colnames(data17.3)[5] <- "P-value"
data17.3$FDR <- NA
data17.3$Species <- "Human"
data17.3$Tissue <- "Whole blood"
data17.3$Disease <- "Never smokers"
data17.3$Article <- "DNA methylation is associated with lung function in never smokers"
data17.3$Journal <- "Respiratory Research"
data17.3$Date <- "2019.12"
data17.3$PMID <- "31791327"
data17.3$Cis_Trans_eqtm <- "Cis"
data17.3$Method <- "Robust linear regression"
data17.3$Sample <- "206"
a <- c(1,2,14,3:13,15,16)
data17.3 <- data17.3[,a]
data17.3 <- data17.3 %>%
  left_join(gtf, by = c("Gene" = "ens"))
a <- c(1,17,3:16)
data17.3 <- data17.3[,a]
colnames(data17.3)[2] <- "GeneSymbol"
has_na <- any(is.na(data17.3$`GeneSymbol`))
na_count <- sum(is.na(data17.3$`GeneSymbol`))
print(has_na)
print(na_count)
data17.3 <- merge(data17.3, 
                  ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                  by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data17.3)[17] <- "CpG_Gene"
colnames(data17.3)[18] <- "Feature"
colnames(data17.3)[19] <- "C_Chr"
colnames(data17.3)[20] <- "C_Pos"
colnames(data17.3)[6] <- "P-value"
data17.3 <- data17.3 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data17.3 <- data17.3[,a]
data17.3 <- data17.3 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(GeneSymbol, CpG, CpG_Gene, CpG_position, everything())
data17.3 <- data17.3[,-c(5,6)]
colnames(data17.3)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data17.3, file = "/data1/renshida/00_final_eqtm//31791327-cis-NTR.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)




data17.4 <- data17[,c(1,3,16,17,18)]
colnames(data17.4)[1] <- "CpG"
colnames(data17.4)[2] <- "Gene"
colnames(data17.4)[3] <- "Beta"
colnames(data17.4)[4] <- "SE"
colnames(data17.4)[5] <- "P-value"
data17.4$FDR <- NA
data17.4$Species <- "Human"
data17.4$Tissue <- "Whole blood"
data17.4$Disease <- "Never smokers"
data17.4$Article <- "DNA methylation is associated with lung function in never smokers"
data17.4$Journal <- "Respiratory Research"
data17.4$Date <- "2019.12"
data17.4$PMID <- "31791327"
data17.4$Cis_Trans_eqtm <- "Cis"
data17.4$Method <- "Robust linear regression"
data17.4$Sample <- "150"
a <- c(1,2,14,3:13,15,16)
data17.4 <- data17.4[,a]
data17.4 <- data17.4 %>%
  left_join(gtf, by = c("Gene" = "ens"))
a <- c(1,17,3:16)
data17.4 <- data17.4[,a]
colnames(data17.4)[2] <- "GeneSymbol"
has_na <- any(is.na(data17.4$`GeneSymbol`))
na_count <- sum(is.na(data17.4$`GeneSymbol`))
print(has_na)
print(na_count)
data17.4 <- merge(data17.4, 
                  ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                  by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data17.4)[17] <- "CpG_Gene"
colnames(data17.4)[18] <- "Feature"
colnames(data17.4)[19] <- "C_Chr"
colnames(data17.4)[20] <- "C_Pos"
colnames(data17.4)[6] <- "P-value"
data17.4 <- data17.4 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data17.4 <- data17.4[,a]
data17.4 <- data17.4 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(GeneSymbol, CpG, CpG_Gene, CpG_position, everything())
data17.4 <- data17.4[,-c(5,6)]
colnames(data17.4)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data17.4, file = "/data1/renshida/00_final_eqtm/31791327-cis-RS.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


#####37598461
data18 <- read.table("37598461-cis.txt",header = T,sep="")
data18 <- data18 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data18 <- data18[,a]
data18 <- data18 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(Gene, CpG, CpG_Gene, CpG_position, everything())
data18 <- data18[,-c(5,6)]
colnames(data18)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
write.table(data18, file = "/data1/renshida/00_final_eqtm/37598461-cis.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)


####33752734
#BiocManager::install("hgu133plus2.db")  # 示例：U133 Plus 2.0 Array
#library(hgu133plus2.db)
data19 <- read_excel("/home/renshida/eqtm/rawdata/human/33752734-cis-trans.xlsx",sheet = "Table S2")
colnames(data19) <- data19[1,]
data19 <- data19[-1,]
data19 <- data19[,c(1,2,5,11)]
colnames(data19)[1] <- "CpG"
colnames(data19)[2] <- "Gene"
colnames(data19)[3] <- "P-value"
data19$Beta <- NA
data19$SE <- NA
data19$FDR <- NA
data19$Species <- "Human"
data19$Tissue <- "Whole blood"
data19$Disease <- "Normal"
data19$Article <- "Epigenome-wide association study of whole blood gene expression in Framingham Heart Study participants provides molecular insight into the potential role of CHRNA5 in cigarette smoking-related lung diseases"
data19$Journal <- "Clinical Epigenetics"
data19$Date <- "2021.03"
data19$PMID <- "33752734"
data19$Cis_Trans_eqtm <- "Cis"
data19$Method <- "Linear mixed model "
data19$Sample <- "4170"
a <- c(1,2,15,3:14,16,17)
data19 <- data19[,a]
a <- c(1,2,3,6,7,4,5,8:17)
data19 <- data19[,a]
a <- c(1:6,8,7,9:17)
data19 <- data19[,a]
str(data19$`P-value`)
data19$`P-value` <- as.numeric(data19$`P-value`)
data19$`P-value` <- 10^(data19$`P-value`)
unique_gene_count <- data19 %>%
  filter(GeneSymbol == "" | is.na(GeneSymbol)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(Gene) %>%                                 # 去重 Gene 列
  nrow() 
unique_gene_count
length(unique(data19$Gene))
a <- c(1,8,3:7,9:17)
data19 <- data19[,a]
data19 <- merge(data19, 
                ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data19)[17] <- "CpG_Gene"
colnames(data19)[18] <- "Feature"
colnames(data19)[19] <- "C_Chr"
colnames(data19)[20] <- "C_Pos"
colnames(data19)[6] <- "P-value"
data19 <- data19 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data19 <- data19[,a]
data19 <- data19 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(GeneSymbol, CpG, CpG_Gene, CpG_position, everything())
data19 <- data19[,-c(5,6)]
colnames(data19)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
data19$`Gene Symbol` <- str_replace(data19$`Gene Symbol`, "\\|.*", "")
data19 <- data19[!is.na(data19$`Gene Symbol`), ]
range(data19$`P-value`)


data19.1 <- read_excel("/home/renshida/eqtm/rawdata/human/33752734-cis-trans.xlsx",sheet = "Table S3")
colnames(data19.1) <- data19.1[1,]
data19.1 <- data19.1[-1,]
data19.1 <- data19.1[,c(1,2,5,11)]
colnames(data19.1)[1] <- "CpG"
colnames(data19.1)[2] <- "Gene"
colnames(data19.1)[3] <- "P-value"
data19.1$Beta <- NA
data19.1$SE <- NA
data19.1$FDR <- NA
data19.1$Species <- "Human"
data19.1$Tissue <- "Whole blood"
data19.1$Disease <- "Normal"
data19.1$Article <- "Epigenome-wide association study of whole blood gene expression in Framingham Heart Study participants provides molecular insight into the potential role of CHRNA5 in cigarette smoking-related lung diseases"
data19.1$Journal <- "Clinical Epigenetics"
data19.1$Date <- "2021.03"
data19.1$PMID <- "33752734"
data19.1$Cis_Trans_eqtm <- "Trans"
data19.1$Method <- "Linear mixed model "
data19.1$Sample <- "4170"
a <- c(1,2,15,3:14,16,17)
data19.1 <- data19.1[,a]
a <- c(1,2,3,6,7,4,5,8:17)
data19.1 <- data19.1[,a]
a <- c(1:6,8,7,9:17)
data19.1 <- data19.1[,a]
str(data19.1$`P-value`)
data19.1$`P-value` <- as.numeric(data19.1$`P-value`)
data19.1$`P-value` <- 10^(data19.1$`P-value`)
unique_gene_count <- data19.1 %>%
  filter(GeneSymbol == "" | is.na(GeneSymbol)) %>%  # 过滤 Genesymbol 为空或 NA 的行
  distinct(Gene) %>%                                 # 去重 Gene 列
  nrow() 
unique_gene_count
length(unique(data19.1$Gene))
a <- c(1,8,3:7,9:17)
data19.1 <- data19.1[,a]
data19.1 <- merge(data19.1, 
                  ann_450k[, c("IlmnID", "UCSC_RefGene_Name", "UCSC_RefGene_Group", "hg38_chr", "hg38_pos")], 
                  by.x = "CpG", by.y = "IlmnID", all.x = TRUE)
colnames(data19.1)[17] <- "CpG_Gene"
colnames(data19.1)[18] <- "Feature"
colnames(data19.1)[19] <- "C_Chr"
colnames(data19.1)[20] <- "C_Pos"
colnames(data19.1)[6] <- "P-value"
data19.1 <- data19.1 %>%
  left_join(ann_450k %>% select(IlmnID, Relation_to_UCSC_CpG_Island), by = c("CpG" = "IlmnID"))
a <- c(2,1,17,19,20,21,18,3:16)
data19.1 <- data19.1[,a]
data19.1 <- data19.1 %>%
  mutate(CpG_position = paste(C_Chr, C_Pos, sep = ":")) %>%
  select(GeneSymbol, CpG, CpG_Gene, CpG_position, everything())
data19.1 <- data19.1[,-c(5,6)]
colnames(data19.1)[1:11] <- c("Gene Symbol", "CpG site","CpG Gene","CpG position","Relation to CpG Island", "RefGene Group","Cis/Trans","Beta","SE","P-value","FDR")
data19.1$`Gene Symbol` <- str_replace(data19.1$`Gene Symbol`, "\\|.*", "")
data19.1 <- data19.1[!is.na(data19.1$`Gene Symbol`), ]
range(data19.1$`P-value`)

data19.2 <- rbind(data19,data19.1)
write.table(data19.2, file = "/data1/renshida/00_final_eqtm//33752734-all.txt", sep = "\t", row.names = FALSE,col.names = T, quote = F)



#####39420233
data20 <- read.csv("//home/renshida/eqtm/rawdata/human/39420233", header = TRUE, sep = ",")
data20 <- data20[,-1]
colnames(data20)[1] <- "CpG"
colnames(data20)[2] <- "Gene"
colnames(data20)[3] <- "P-value"
colnames(data20)[5] <- "Beta"
data20 <- data20[,c(1,2,5,3,4)]
data20$Species <- "Human"
data20$Tissue <- "Placenta"
data20$Disease <- "Imprint control regions"
data20$Article <- "In-depth characterization of the placental imprintome reveals novel differentially methylated regions across birth weight categories"
data20$Journal <- "Epigenetics"
data20$Date <- "2020.01-02"
data20$PMID <- "31403346"
data20$Cis_Trans_eqtm <- "Cis"
data20$Method <- "Linear regression model"
data20$Sample <- "163"
data20$SE <- NA
a <- c(1:2, 13, 3:12, 14:ncol(data20))
data20 <- data20[, a]
a <- c(1:4, 16, 5:15)
data20 <- data20[, a]
write.table(data20, file = "/home/renshida/eqtm/data/human/31403346-cis.txt", sep = "\t", row.names = FALSE)










































