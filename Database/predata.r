setwd("/home/renshida/eqtm/rawdata/human")
getwd()
#install.packages("dplyr")
#install.packages("readxl")
#install.packages("readr")
#install.packages("tidyr")
# # BiocManager::install("biomaRt")
# BiocManager::install("dbplyr")
# BiocManager::install("BiocFileCache")
# library(BiocManager)
# library(biomaRt)
#library(BiocLite)
if (!requireNamespace("AnnotationDbi", quietly = TRUE)) {
  install.packages("AnnotationDbi")
}
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  install.packages("org.Hs.eg.db")
}
library(AnnotationDbi)
library(org.Hs.eg.db)
library(tidyr)
library(readr)
library(readxl)
library(dplyr)



#######33092652
file_path <- "33092652-cis-trans.xlsx"
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
combined_data <- bind_rows(processed_data_list)
colnames(combined_data)[1] <- "CpG"
colnames(combined_data)[2] <- "Beta"
colnames(combined_data)[3] <- "P-value"
a <- c(1,8,2:7,9:ncol(combined_data))
combined_data <- combined_data[, a]
combined_data <- combined_data[,c(1,2,3,4)]
combined_data$SE <- NA
combined_data$FDR <- NA
combined_data$Species <- "Human"
combined_data$Tissue <- "Whole blood"
combined_data$Disease <- "Differences between smokers and never smokers"
combined_data$Article <- "Smoking-related changes in DNA methylation and gene expression are associated with cardio-metabolic traits"
combined_data$Journal <- "Clinical Epigenetics"
combined_data$Date <- "2020.10"
combined_data$PMID <- "33092652"
combined_data$Cis_Trans_eqtm <- NA
combined_data$Method <- "Linear mixed model "
combined_data$Sample <- "716"
#na_columns <- apply(combined_data, 2, function(x) any(is.na(x)))
#print(na_columns)
a <- c(1:3,5,4,6:ncol(combined_data))
combined_data <- combined_data[, a]
combined_data$'P-value' <- as.numeric(combined_data$`P-value`)
pvalue_range <- range(combined_data$'P-value', na.rm = TRUE)
print(pvalue_range)
a <- c(1:2,14,3:13,15,16)
combined_data <- combined_data[,a]
#a <- c(1:4,6,5,7:ncol(combined_data))
#combined_data <- combined_data[, a]
#colnames(combined_data)[5] <- "P-value"
#colnames(combined_data)[6] <- "FDR"
write.table(combined_data, file = "/home/renshida/eqtm/data/human/33092652-cis-trans-rot.txt", sep = "\t", row.names = FALSE)



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
#a <- c(1:4,6,5,7:ncol(combined_data1))
#combined_data1 <- combined_data1[, a]
#colnames(combined_data1)[5] <- "P-value"
#colnames(combined_data1)[6] <- "FDR"
a <- c(1:2,14,3:13,15,16)
combined_data1 <- combined_data1[,a]
write.table(combined_data1, file = "/home/renshida/eqtm/data/human/33092652-cis-trans-kora.txt", sep = "\t", row.names = FALSE)
range(combined_data1$FDR)

################37271320
rawdata <- read_excel("37271320-cis.xlsx")
rawdata <- rawdata[, c(1, 2, 6, 4, 5)]
colnames(rawdata)[1] <- "CpG"
colnames(rawdata)[2] <- "Gene"
colnames(rawdata)[3] <- "Beta"
colnames(rawdata)[4] <- "P-value"
colnames(rawdata)[5] <- "FDR"
rawdata$Species <- "Human"
rawdata$Tissue <- "Nasal epithelial"
rawdata$Disease <- "Atopic asthma"
rawdata$Article <- "Cis- and trans-eQTM analysis reveals novel epigenetic and transcriptomic immune markers of atopic asthma in airway epithelium"
rawdata$Journal <- "Journal Of Allergy And Clinical Immunology"
rawdata$Date <- "2023.10"
rawdata$PMID <- "37271320"
rawdata$Method <- "Linear regression model"
rawdata$Cis_Trans_eqtm <- "Cis"
rawdata$Sample <- "258"
a <- c(1:2, 14, 3:13, 15:ncol(rawdata))
rawdata <- rawdata[, a]
rawdata$SE <- NA
a <- c(1:4,16,5:15)
rawdata <- rawdata[, a]
write.table(rawdata, file = "/home/renshida/eqtm/data/human/37271320-cis.txt", sep = "\t", row.names = FALSE)






