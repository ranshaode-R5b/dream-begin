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










































