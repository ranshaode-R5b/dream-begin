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


rawdata <- read_excel("37271320-trans1.xlsx")
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
rawdata$Cis_Trans_eqtm <- "Trans"
rawdata$Sample <- "258"
a <- c(1:2, 14, 3:13, 15:ncol(rawdata))
rawdata <- rawdata[, a]
rawdata$SE <- NA
a <- c(1:4,16,5:15)
rawdata <- rawdata[, a]
write.table(rawdata, file = "/home/renshida/eqtm/data/human/37271320-trans1.txt", sep = "\t", row.names = FALSE)



rawdata <- read_excel("37271320-trans2.xlsx")
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
rawdata$Cis_Trans_eqtm <- "Trans"
rawdata$Sample <- "258"
a <- c(1:2, 14, 3:13, 15:ncol(rawdata))
rawdata <- rawdata[, a]
rawdata$SE <- NA
a <- c(1:4,16,5:15)
rawdata <- rawdata[, a]
write.table(rawdata, file = "/home/renshida/eqtm/data/human/37271320-trans2.txt", sep = "\t", row.names = FALSE)



rawdata <- read_excel("37271320-trans3.xlsx")
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
rawdata$Cis_Trans_eqtm <- "Trans"
rawdata$Sample <- "258"
a <- c(1:2, 14, 3:13, 15:ncol(rawdata))
rawdata <- rawdata[, a]
rawdata$SE <- NA
a <- c(1:4,16,5:15)
rawdata <- rawdata[, a]
write.table(rawdata, file = "/home/renshida/eqtm/data/human/37271320-trans3.txt", sep = "\t", row.names = FALSE)
range(rawdata$pvalue)




###########35426765
data_E7 <- read_excel("35426765-cis.xlsx", sheet = "Table E7", col_names = FALSE)
colnames(data_E7) <- data_E7[2, ]
data_E7 <- data_E7[-1, ]
data_E7 <- data_E7[, c(1,7,2,3,5,8)]
colnames(data_E7)[1] <- "CpG"
colnames(data_E7)[2] <- "Gene"
colnames(data_E7)[5] <- "P-value"
colnames(data_E7)[6] <- "FDR"
data_E7$Species <- "Human"
data_E7$Tissue <- "Nasal epithelial"
data_E7$Disease <- "HIV and chronic obstructive pulmonary disease"
data_E7$Article <- "Airway Aging and Methylation Disruptions in HIV-associated Chronic Obstructive Pulmonary Disease"
data_E7$Journal <- "American Journal of Respiratory and Critical Care Medicine"
data_E7$Date <- "2022.07"
data_E7$PMID <- "35426765"
data_E7$Method <- "Robust linear model"
data_E7$Cis_Trans_eqtm <- "Cis"
data_E7$Sample <- "76"
a <- c(1:2, 15, 3:14, 16:ncol(data_E7))
data_E7 <- data_E7[, a]
data_E7$FDR <- as.numeric(data_E7$FDR)
pvalue_range <- range(data_E7$FDR, na.rm = TRUE)
print(pvalue_range)
write.table(data_E7, file = "/home/renshida/eqtm/data/human/35426765-cis.txt", sep = "\t", row.names = FALSE)



#######31403346
data1 <- read.csv("31403346-cis.csv", header = TRUE, sep = ",")
data1 <- data1[,-1]
colnames(data1)[1] <- "CpG"
colnames(data1)[2] <- "Gene"
colnames(data1)[3] <- "P-value"
colnames(data1)[5] <- "Beta"
data1 <- data1[,c(1,2,5,3,4)]
data1$Species <- "Human"
data1$Tissue <- "Placenta"
data1$Disease <- "Imprint control regions"
data1$Article <- "In-depth characterization of the placental imprintome reveals novel differentially methylated regions across birth weight categories"
data1$Journal <- "Epigenetics"
data1$Date <- "2020.01-02"
data1$PMID <- "31403346"
data1$Cis_Trans_eqtm <- "Cis"
data1$Method <- "Linear regression model"
data1$Sample <- "163"
data1$SE <- NA
a <- c(1:2, 13, 3:12, 14:ncol(data1))
data1 <- data1[, a]
a <- c(1:4, 16, 5:15)
data1 <- data1[, a]
write.table(data1, file = "/home/renshida/eqtm/data/human/31403346-cis.txt", sep = "\t", row.names = FALSE)




######32493484
data2 <- read_excel("32493484-cis.xlsx",sheet = "Sheet2")
colnames(data2)[1] <- "CpG"
colnames(data2)[2] <- "Gene"
colnames(data2)[3] <- "Beta"
colnames(data2)[5] <- "P-value"
colnames(data2)[6] <- "FDR"
data2$Species <- "Human"
data2$Tissue <- "Placenta"
data2$Disease <- "Birthweight"
data2$Article <- "DNA methylation loci in placenta associated with birthweight and expression of genes relevant for early development and adult diseases "
data2$Journal <- "Clinical Epigenetics"
data2$Date <- "2020.06"
data2$PMID <- "32493484"
data2$Cis_Trans_eqtm <- "Cis"
data2$Method <- "Linear regression model"
data2$Sample <- "75"
a <- c(1:2, 14, 3:13, 15:ncol(data2))
data2 <- data2[, a]
write.table(data2, file = "/home/renshida/eqtm/data/human/32493484-cis.txt", sep = "\t", row.names = FALSE)



#######
data3 <- read_excel("/home/renshida/eqtm/rawdata/human/37598461-cis.xlsx",col_names = FALSE)
data3 <- data3[-1,]
colnames(data3) <- data3[1,]
data3 <- data3[-1,]
data3 <- data3[,c(1,8,9,10)]
colnames(data3)[1] <- "CpG"
colnames(data3)[2] <- "Gene"
colnames(data3)[3] <- "Beta"
colnames(data3)[4] <- "P-value"
data3$Species <- "Human"
data3$SE <- NA
data3$FDR <- NA
data3$Tissue <- "Whole blood"
data3$Disease <- "Serum immunoglobulin E"
data3$Article <- "Epigenome-wide DNA methylation association study of circulating IgE levels identifies novel targets for asthma"
data3$Journal <- "EBioMedicine"
data3$Date <- "2023.09"
data3$PMID <- "37598461"
data3$Cis_Trans_eqtm <- "Cis"
data3$Method <- "Linear regression model"
data3$Sample <- "1045"
a <- c(1:2,14,3:13,15:ncol(data3))
data3 <- data3[, a]
a <- c(1:4,7,5:6,8:ncol(data3))
data3 <- data3[, a]
a <- c(1:6,8,7,9:ncol(data3))
data3 <- data3[, a]
#a <- c(1:5,7,6,8:ncol(data3))
#data3 <- data3[, a]
#colnames(data3)[6] <- "P-value"
#colnames(data3)[7] <- "FDR"
write.table(data3, file = "/home/renshida/eqtm/data/human/37598461-cis.txt", sep = "\t", row.names = FALSE)
data3$FDR <- as.numeric(data3$FDR)
pvalue_range <- range(data3$FDR, na.rm = TRUE)
print(pvalue_range)






#######
data4 <- read_excel("31282290-cis-mir.xlsx",col_names = FALSE)
data4 <- data4[-1,]
colnames(data4) <- data4[1,]
data4 <- data4[-1,]
data4 <- data4[,c(1,11,3,6,7)]
colnames(data4)[1] <- "CpG"
colnames(data4)[2] <- "Gene"
colnames(data4)[3] <- "Beta"
colnames(data4)[4] <- "P-value"
data4$Species <- "Human"
data4$SE <- NA
data4$Tissue <- "Whole blood"
data4$Disease <- "MicroRNAs"
data4$Article <- "Epigenome-wide association study of DNA methylation and microRNA expression highlights novel pathways for human complex traits"
data4$Journal <- "Epigenetics"
data4$Date <- "2020.01-02"
data4$PMID <- "31282290"
data4$Cis_Trans_eqtm <- "Cis"
data4$Method <- "Linear mixed model"
data4$Sample <- "3345"
a <- c(1:3,7,4:6,8:ncol(data4))
data4 <- data4[, a]
data4$'FDR' <- as.numeric(data4$`FDR`)
pvalue_range <- range(data4$'FDR', na.rm = TRUE)
print(pvalue_range)
a <- c(1,2,14,3:13,15,16)
data4 <- data4[,a]
write.table(data4, file = "/home/renshida/eqtm/data/human/31282290-cis-mir.txt", sep = "\t", row.names = FALSE)




######
data5 <- read_excel("37047575-cis-region.xlsx")
data5 <- data5 %>% unite(DMR, c("seqnames", "start", "end"), sep = "_")
data5 <- data5[,c(1,9,4,6,10)]
colnames(data5)[1] <- "DMR"
colnames(data5)[3] <- "Gene"
colnames(data5)[4] <- "Beta"
colnames(data5)[5] <- "P-value"
data5$Species <- "Human"
data5$SE <- NA
data5$FDR <- NA
data5$Tissue <- "Whole blood"
data5$Disease <- "Fetal alcohol spectrum disorder"
data5$Article <- "Expression Quantitative Trait Methylation Analysis Identifies Whole Blood Molecular Footprint in Fetal Alcohol Spectrum Disorder (FASD)"
data5$Journal <- "International Journal of Molecular Sciences"
data5$Date <- "2023.04"
data5$PMID <- "37047575"
data5$Cis_Trans_eqtm <- "Cis"
data5$Method <- "Pearson’s correlation coefficient"
data5$Sample <- "63"
a <- c(1:3,15,4:14,16,17)
data5 <- data5[,a]
a <- c(1:5,8,6:7,9:ncol(data5))
data5 <- data5[,a]
a <- c(1:7,9,8,10:17)
data5 <- data5[,a]
write.table(data5, file = "/home/renshida/eqtm/data/human/37047575-cis-DMR.txt", sep = "\t", row.names = FALSE)



data6 <- read_tsv("29914364-GTP.txt") 
range(data6$p.val)
data6$p.val <- as.numeric(data6$p.val)
pvalue_range <- range(data6$p.val, na.rm = TRUE)
print(pvalue_range)
range(data6$beta)
data6 <- data6[,c(2,7,15,13)]
colnames(data6)[1] <- "CpG"
colnames(data6)[2] <- "Gene"
colnames(data6)[3] <- "Beta"
colnames(data6)[4] <- "P-value"
data6$FDR <- NA
data6$Species <- "Human"
data6$SE <- NA
data6$Tissue <- "Whole blood"
data6$Disease <- "Post traumatic stress disorder"
data6$Article <- "An integrated -omics analysis of the epigenetic landscape of gene expression in human blood cells"
data6$Journal <- "BMC Genomics"
data6$Date <- "2018.06"
data6$PMID <- "29914364"
data6$Cis_Trans_eqtm <- NA
data6$Method <- "Linear mixed model"
data6$Sample <- "333"
a <- c(1:3,7,4:6,8:ncol(data6))
data6 <- data6[,a]
a <- c(1:2,14,3:13,15:ncol(data6))
data6 <- data6[,a]
write.table(data6, file = "/home/renshida/eqtm/data/human/29914364-GTP.txt", sep = "\t", row.names = FALSE)



data7 <- read_tsv("29914364-MESA.txt") 
range(data7$p.val)
data7$p.val <- as.numeric(data7$p.val)
pvalue_range <- range(data7$p.val, na.rm = TRUE)
print(pvalue_range)
range(data7$beta)
data7 <- data7[,c(1,7,15,13)]
colnames(data7)[1] <- "CpG"
colnames(data7)[2] <- "Gene"
colnames(data7)[3] <- "Beta"
colnames(data7)[4] <- "P-value"
data7$FDR <- NA
data7$Species <- "Human"
data7$SE <- NA
data7$Tissue <- "Monocytes"
data7$Disease <- "Atherosclerosis"
data7$Article <- "An integrated -omics analysis of the epigenetic landscape of gene expression in human blood cells"
data7$Journal <- "BMC Genomics"
data7$Date <- "2018.06"
data7$PMID <- "29914364"
data7$Cis_Trans_eqtm <- NA
data7$Method <- "Linear mixed model"
data7$Sample <- "1202"
a <- c(1:3,7,4:6,8:ncol(data7))
data7 <- data7[,a]
a <- c(1:2,14,3:13,15:ncol(data6))
data7 <- data7[,a]
write.table(data7, file = "/home/renshida/eqtm/data/human/29914364-MESA.txt", sep = "\t", row.names = FALSE)
range(data7$`P-value`)





# data8 <- read_excel("39009577-cis.xlsx")
# data8 <- data8[,c(2,17,8,20,1,22)]
# data8 <- data8 %>%
#   separate(`Beta (SE)`, into = c("Beta1", "Beta2", "Beta3", "Beta4"), sep = ";\\s*", convert = TRUE)
# data8 <- data8 %>%
#   separate(Beta1, into = c("Beta1_value", "Beta1_SE"), sep = "\\(", convert = TRUE) %>%
#   separate(Beta2, into = c("Beta2_value", "Beta2_SE"), sep = "\\(", convert = TRUE) %>%
#   separate(Beta3, into = c("Beta3_value", "Beta3_SE"), sep = "\\(", convert = TRUE) %>%
#   separate(Beta4, into = c("Beta4_value", "Beta4_SE"), sep = "\\(", convert = TRUE)
# data8 <- data8 %>%
#   mutate(across(everything(), ~ gsub("\\)", "", .)))
# colnames(data8)[1] <- "CpG"
# colnames(data8)[2] <- "Gene"
# colnames(data8)[3] <- "Cis_Trans_eqtm"
# data8RS <- data8[,c(1,2,3,4,5,12,13)]
# colnames(data8RS)[4] <- "Beta"
# colnames(data8RS)[5] <- "SE"
# colnames(data8RS)[6] <- "P-value"
# data8RS$Species <- "Human"
# data8RS$Tissue <- "whole blood"
# data8RS$Disease <- "Depression"
# data8RS$Article <- "Decoding depression: a comprehensive multi-cohort exploration of blood DNA methylation using machine learning and deep learning approaches"
# data8RS$Journal <- "Translational Psychiatry"
# data8RS$Date <- "2024.07"
# data8RS$PMID <- "39009577"
# data8RS$Method <- "Linear mixed model"
# data8RS$Sample <- "1202"





########
data9c <- read_excel("37233989-cis.xlsx")
colnames(data9c) <- data9c[1,]
data9c <- data9c[-1,]
data9c <- data9c[,c(1,9,11,13,14)]
colnames(data9c)[1] <- "CpG"
colnames(data9c)[2] <- "Gene"
colnames(data9c)[3] <- "Beta"
colnames(data9c)[4] <- "P-value"
colnames(data9c)[5] <- "FDR"
data9c$`P-value` <- 10^ as.numeric(data9c$`P-value`)
data9c$`FDR` <- 10^ as.numeric(data9c$`FDR`)
data9c$Species <- "Human"
data9c$SE <- NA
data9c$Tissue <- "Whole blood"
data9c$Disease <- "Vitamin C intake"
data9c$Article <- "Dietary and supplemental intake of vitamins C and E is associated with altered DNA methylation in an epigenome-wide association study meta-analysis"
data9c$Journal <- "Epigenetics"
data9c$Date <- "2023.12"
data9c$PMID <- "37233989"
data9c$Cis_Trans_eqtm <- "Cis"
data9c$Method <- "Linear mixed model"
data9c$Sample <- "3860"
a <- c(1:3,7,4:6,8:ncol(data9c))
data9c <- data9c[,a]
a <- c(1:2,14,3:13,15:ncol(data9c))
data9c <- data9c[,a]
write.table(data9c, file = "/home/renshida/eqtm/data/human/37233989-cis-vc.txt", sep = "\t", row.names = FALSE)


data9e <- read_excel("37233989-cis.xlsx",sheet = "ST11")
colnames(data9e) <- data9e[1,]
data9e <- data9e[-1,]
data9e <- data9e[,c(1,9,11,13,14)]
colnames(data9e)[1] <- "CpG"
colnames(data9e)[2] <- "Gene"
colnames(data9e)[3] <- "Beta"
colnames(data9e)[4] <- "P-value"
colnames(data9e)[5] <- "FDR"
data9e$`P-value` <- 10^ as.numeric(data9e$`P-value`)
data9e$`FDR` <- 10^ as.numeric(data9e$`FDR`)
data9e$Species <- "Human"
data9e$SE <- NA
data9e$Tissue <- "Whole blood"
data9e$Disease <- "Vitamin E intake"
data9e$Article <- "Dietary and supplemental intake of vitamins C and E is associated with altered DNA methylation in an epigenome-wide association study meta-analysis"
data9e$Journal <- "Epigenetics"
data9e$Date <- "2023.12"
data9e$PMID <- "37233989"
data9e$Cis_Trans_eqtm <- "Cis"
data9e$Method <- "Linear mixed model"
data9e$Sample <- "3860"
a <- c(1:3,7,4:6,8:ncol(data9e))
data9e <- data9e[,a]
a <- c(1:2,14,3:13,15:ncol(data9e))
data9e <- data9e[,a]
write.table(data9e, file = "/home/renshida/eqtm/data/human/37233989-cis-ve.txt", sep = "\t", row.names = FALSE)




########
# data10 <- read_tsv("/home/renshida/eqtm/rawdata/human/35302492-cis-adj.txt")
# head(data10)
# colnames(data10)[1] <- "CpG"
# colnames(data10)[2] <- "Gene"
# colnames(data10)[3] <- "Beta"
# colnames(data10)[4] <- "SE"
# colnames(data10)[5] <- "P-value"
# data10 <- data10[,c(1,2,3,4,5)]
# data10$Species <- "Human"
# data10$Tissue <- "Whole blood"
# data10$Disease <- "Autosomal"
# data10$Article <- "Identification of autosomal cis expression quantitative trait methylation (cis eQTMs) in children's blood"
# data10$Journal <- "Elife"
# data10$Date <- "2022.03"
# data10$PMID <- "35302492"
# data10$Cis_Trans_eqtm <- "Cis"
# data10$Method <- "Linear regression model"
# data10$Sample <- "832"
# data10$FDR <- NA
# data10$'P-value' <- as.numeric(data10$'P-value')
# pvalue_range <- range(data10$'P-value', na.rm = TRUE)
# print(pvalue_range)
# a <- c(1:2,13,3:12,14:16)
# data10 <- data10[,a]
# a <- c(1:6,16,7:15)
# data10 <- data10[,a]
# xxx <- read_tsv("/home/renshida/eqtm/rawdata/human/35302492-cis-adj.txt")
# data10$Genesymbol <- xxx$TC_gene
# xxx <- NULL
# write.table(data10, file = "/home/renshida/eqtm/data/human/35302492-cis-adj.txt", sep = "\t", row.names = FALSE)
# write.table(data10[,c(2,17)], file = "/home/renshida/eqtm/data/trans_g_id/35302492-cis-adj.txt", sep = "\t", row.names = FALSE)
# has_na <- any(is.na(data10$TC_gene))
# na_count <- sum(is.na(data10$TC_gene))
# print(has_na)
# print(na_count)



data11 <- read_tsv("35302492-cis-unadj.txt")
colnames(data11)[1] <- "CpG"
colnames(data11)[3] <- "Beta"
colnames(data11)[4] <- "SE"
colnames(data11)[5] <- "P-value"
data11 <- data11[,c(1,15,3,4,5)]
colnames(data11)[2] <- "Gene"
data11$Species <- "Human"
data11$Tissue <- "Whole blood"
data11$Disease <- "Autosomal"
data11$Article <- "Identification of autosomal cis expression quantitative trait methylation (cis eQTMs) in children's blood"
data11$Journal <- "Elife"
data11$Date <- "2022.03"
data11$PMID <- "35302492"
data11$Cis_Trans_eqtm <- "Cis"
data11$Method <- "Linear regression model"
data11$Sample <- "832"
data11$FDR <- NA
data11$p.value <- as.numeric(data11$p.value)
a <- c(1:2,13,3:12,14:16)
data11 <- data11[,a]
a <- c(1:6,16,7:15)
data11 <- data11[,a]
data11 <- data11[!is.na(data11$Gene), ]
write.table(data11, file = "/home/renshida/eqtm/data/human/35302492-cis-unadj.txt", sep = "\t", row.names = FALSE)
has_na <- any(is.na(data11$Gene))
na_count <- sum(is.na(data11$Gene))
print(has_na)
print(na_count)




#####
data12 <- read_excel("39393618-cis.xlsx")
data12 <-  data12[,c(1,10,3,5,6)]
colnames(data12)[1] <- "CpG"
colnames(data12)[2] <- "Gene"
colnames(data12)[3] <- "Beta"
colnames(data12)[4] <- "P-value"
colnames(data12)[5] <- "FDR"
data12$Species <- "Human"
data12$Tissue <- "Peripheral blood"
data12$Disease <- "Glucocorticoid receptor activation"
data12$Article <- "Multiomics Analysis of the Molecular Response to Glucocorticoids: Insights Into Shared Genetic Risk From Psychiatric to Medical Disorders"
data12$Journal <- "Biological Psychiatry"
data12$Date <- "2024.10 "
data12$PMID <- "39393618"
data12$Cis_Trans_eqtm <- "Cis"
data12$Method <- "Linear regression model"
data12$Sample <- "398"
data12$SE <- NA
a <- c(1:2,13,3:12,14:16)
data12 <- data12[,a]
a <- c(1:4,16,5:15)
data12 <- data12[,a]
write.table(data12, file = "/home/renshida/eqtm/data/human/39393618-cis.txt", sep = "\t", row.names = FALSE)




#########
data13 <- read_tsv("38438351-cis.txt")
data13 <- data13[,c(2,1,4,5,3,9)]
colnames(data13)[1] <- "CpG"
colnames(data13)[2] <- "Gene"
colnames(data13)[3] <- "Beta"
colnames(data13)[4] <- "SE"
colnames(data13)[5] <- "P-value"
colnames(data13)[6] <- "FDR"
data13$Species <- "Human"
data13$Tissue <- "Retina"
data13$Disease <- "Normal"
data13$Article <- "QTL mapping of human retina DNA methylation identifies 87 gene-epigenome interactions in age-related macular degeneration"
data13$Journal <- "Nature Communications"
data13$Date <- "2024.03"
data13$PMID <- "38438351"
data13$Cis_Trans_eqtm <- "Cis"
data13$Method <- "Linear regression model"
data13$Sample <- "152"
a <- c(1:2,14,3:13,15,16)
data13 <- data13[,a]
data13$`P-value` <- as.numeric(data13$`P-value`)
data13$FDR <- as.numeric(data13$FDR)
range(data13$`P-value`)
range(data13$FDR)
##转cpg.r
write.table(data13, file = "/home/renshida/eqtm/data/human/38438351-cis.txt", sep = "\t", row.names = FALSE)





##########
data14 <- read_excel("32928872-cis-cite-region.xlsx")
colnames(data14) <- data14[5,]
data14 <- data14[-c(1:5),]
data14 <- data14[,c(1,6,7,8,9,10)]
colnames(data14)[1] <- "CpG"
colnames(data14)[2] <- "Gene"
colnames(data14)[3] <- "Beta"
colnames(data14)[4] <- "SE"
colnames(data14)[5] <- "P-value"
colnames(data14)[6] <- "FDR"
data14$Species <- "Human"
data14$Tissue <- "Adipose"
data14$Disease <- "Normal"
data14$Article <- "Integrative Analysis of Glucometabolic Traits, Adipose Tissue DNA Methylation, and Gene Expression Identifies Epigenetic Regulatory Mechanisms of Insulin Resistance and Obesity in African Americans"
data14$Journal <- "Diabetes"
data14$Date <- "2020.12"
data14$PMID <- "32928872"
data14$Cis_Trans_eqtm <- "Cis"
data14$Method <- "Linear regression model"
data14$Sample <- "230"
a <- c(1:2,14,3:13,15,16)
data14 <- data14[,a]
write.table(data14, file = "/home/renshida/eqtm/data/human/32928872-cis-cite.txt", sep = "\t", row.names = FALSE)


data14a <- read_excel("32928872-cis-cite-region.xlsx",sheet = "table-S2B")
colnames(data14a) <- data14a[6,]
data14a <- data14a[-c(1:6),]
data14a <- data14a[,c(1,22,6,7,8,9,10)]
colnames(data14a)[1] <- "CpG region"
colnames(data14a)[2] <- "nCpG"
colnames(data14a)[3] <- "Gene"
colnames(data14a)[4] <- "Beta"
colnames(data14a)[5] <- "SE"
colnames(data14a)[6] <- "P-value"
colnames(data14a)[7] <- "FDR"
data14a$Species <- "Human"
data14a$Tissue <- "Adipose"
data14a$Disease <- "Normal"
data14a$Article <- "Integrative Analysis of Glucometabolic Traits, Adipose Tissue DNA Methylation, and Gene Expression Identifies Epigenetic Regulatory Mechanisms of Insulin Resistance and Obesity in African Americans"
data14a$Journal <- "Diabetes"
data14a$Date <- "2020.12"
data14a$PMID <- "32928872"
data14a$Cis_Trans_eqtm <- "Cis"
data14a$Method <- "Linear regression model"
data14a$Sample <- "230"
a <- c(1:3,15,4:14,16,17)
data14a <- data14a[,a]
write.table(data14a, file = "/home/renshida/eqtm/data/human/32928872-cis-region.txt", sep = "\t", row.names = FALSE)




#######
data15 <- read_tsv("31076557-cis.txt")
data15 <- data15[!is.na(data15$pv), ]
range(data15$pv)
range(data15$qv_bh)
data15 <- data15[,c(6,3,7,9,10)]
data15 <- data15[!is.na(data15$pv), ]
colnames(data15)[1] <- "CpG"
colnames(data15)[2] <- "Gene"
colnames(data15)[3] <- "Beta"
colnames(data15)[4] <- "P-value"
colnames(data15)[5] <- "FDR"
data15 <- data15 %>% filter(FDR < 0.01)
data15$Species <- "Human"
data15$Tissue <- "Skeletal muscle"
data15$Disease <- "Normal"
data15$Article <- "Integrative analysis of gene expression, DNA methylation, physiological traits, and genetic variation in human skeletal muscle"
data15$Journal <- "Proceedings of the National Academy of Sciences of the United States of America"
data15$Date <- "2019.05"
data15$PMID <- "31076557"
data15$Cis_Trans_eqtm <- "Cis"
data15$Method <- "Linear regression model"
data15$Sample <- 318
data15$SE <- NA
a <- c(1,2,13,3:12,14,15,16)
data15 <- data15[,a]
a <- c(1:4,16,5:15)
data15 <- data15[,a]
write.table(data15, file = "/home/renshida/eqtm/data/human/31076557-cis.txt", sep = "\t", row.names = FALSE)
# library(biomaRt)
# 
# # 使用 Ensembl 数据库
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
# ensembl <- useEnsembl(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", mirror = "asia")
# ensembl <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl", host = "https://useast.ensembl.org")
# # 去除版本号（小数点后面的部分）
# ensembl_ids <- sub("\\..*", "", data15$Gene)
# 
# # 获取基因符号
# gene_info <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
#                    filters = "ensembl_gene_id",
#                    values = ensembl_ids,
#                    mart = ensembl)
# 
# # 将基因符号添加到原始数据框
# data15 <- merge(data15, gene_info, by.x = "Gene", by.y = "ensembl_gene_id", all.x = TRUE)

library(AnnotationDbi)
library(org.Hs.eg.db)
data15$Gene <- sub("\\..*", "", data15$Gene)
data15$GeneSymbol <- mapIds(org.Hs.eg.db,
                            keys = data15$Gene,
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
duplicated_genes <- mapping[duplicated(mapping$ENSEMBL) | duplicated(mapping$ENSEMBL, fromLast = TRUE), ]

# 查看具有多重映射的基因
print(duplicated_genes)

# 处理多对一映射，例如，仅保留第一个基因符号
unique_mapping <- mapping[!duplicated(mapping$ENSEMBL), ]






######
#(data16 <- read_tsv("37980427.txt")
#data16 <- data16[,c(1,3,5,6,7,8)]
# colnames(data16)[1] <- "CpG"
# colnames(data16)[2] <- "Gene"
# colnames(data16)[3] <- "Beta"
# colnames(data16)[4] <- "SE"
# colnames(data16)[5] <- "P-value"
# colnames(data16)[6] <- "FDR"
# data16$Species <- "Human"
# data16$Tissue <- "Kidney"
# data16$Disease <- "Kidney stone"
# data16$Article <- "Integrative genome-wide analyses identify novel loci associated with kidney stones and provide insights into its genetic architecture"
# data16$Journal <- "Nature Communications"
# data16$Date <- "2023.11"
# data16$PMID <- "37980427"
# data16$Cis_Trans_eqtm <- NA
# data16$Method <- "Linear regression model"
# data16$Sample <- 414
# a <- c(1:2,14,3:13,15,16)
# # data16 <- data16[,a]
# write.table(data16, file = "/home/renshida/eqtm/data/human/37980427.txt", sep = "\t", row.names = FALSE)
# data16 <- NULL




############
# data17 <- read_excel("35176128-cis.xlsx")
# colnames(data17) <- data17[1,]
# data17 <- data17[-1,]
# data17p <- data17[,c(2,1,3)]
# colnames(data17p)[1] <- "CpG"
# colnames(data17p)[2] <- "Gene"
# colnames(data17p)[3] <- "Beta"
# data17p <- data17p %>%
#   filter(rowSums(across(everything(), ~ . == "NA")) == 0)
# data17p$Cis_Trans_eqtm <- "Cis"
# data17p$SE <- NA
# data17p$'P-value' <- NA
# data17p$FDR <- NA
# data17p$Species <- "Human"
# data17p$Tissue <- "Skin cutaneous melanoma"
# data17p$Disease <- "Skin cutaneous melanoma"
# data17p$Article <- "Tumor expression quantitative trait methylation screening reveals distinct CpG panels for deconvolving cancer immune signatures"
# data17p$Journal <- "Cancer Research"
# data17p$Date <- "2022.05"
# data17p$PMID <- "35176128"
# data17p$Method <- "Spearman correlation coefficient"
# data17p$Sample <- "471"
# a <- c(1:2,4,3,5:16)
# data17p <- data17p[,a]



# data17n <- data17[,c(5,1,6)]
# colnames(data17n)[1] <- "CpG"
# colnames(data17n)[2] <- "Gene"
# colnames(data17n)[3] <- "Beta"
# data17n <- data17n %>%
#   filter(rowSums(across(everything(), ~ . == "NA")) == 0)
# data17n$Cis_Trans_eqtm <- "Cis"
# data17n$SE <- NA
# data17n$'P-value' <- NA
# data17n$FDR <- NA
# data17n$Species <- "Human"
# data17n$Tissue <- "Skin cutaneous melanoma"
# data17n$Disease <- "Skin cutaneous melanoma"
# data17n$Article <- "Tumor expression quantitative trait methylation screening reveals distinct CpG panels for deconvolving cancer immune signatures"
# data17n$Journal <- "Cancer Research"
# data17n$Date <- "2022.05"
# data17n$PMID <- "35176128"
# data17n$Method <- "Spearman correlation coefficient"
# data17n$Sample <- "471"
# a <- c(1:2,4,3,5:16)
# data17n <- data17n[,a]
# data17a <- rbind(data17n, data17p)
# write.table(data17a, file = "/home/renshida/eqtm/data/human/35176128-cis.txt", sep = "\t", row.names = FALSE)




data18 <- read_excel("30390659-cis.xlsx")
colnames(data18) <- data18[1,]
data18 <- data18[-1,]
data18 <- data18[,c(1,11,2,3,4,5)]
colnames(data18)[1] <- "CpG"
colnames(data18)[2] <- "Gene"
colnames(data18)[3] <- "Beta"
colnames(data18)[4] <- "SE"
colnames(data18)[5] <- "P-value"
colnames(data18)[6] <- "FDR"
data18$Species <- "Human"
data18$Tissue <- "Lung"
data18$Disease <- "Differences between smokers and never smokers"
data18$Article <- "From blood to lung tissue: effect of cigarette smoke on DNA methylation and lung function"
data18$Journal <- "Respiratory Research"
data18$Date <- "2018.11"
data18$PMID <- "30390659"
data18$Cis_Trans_eqtm <- "Cis"
data18$Method <- "Linear regression model"
data18$Sample <- 36
a <- c(1:2,14,3:13,15:16)
data18 <- data18[,a]
data18 <- data18[!is.na(data18$Gene), ]
# xxx <- read_excel("30390659-cis.xlsx")
# colnames(xxx) <- xxx[1,]
# xxx <- xxx[-1,]
# data18$Genesymbol <- xxx$Genename
# xxx <- NULL
# write.table(data18[,c(2,17)], file = "/home/renshida/eqtm/data/trans_g_id/30390659-cis.txt", sep = "\t", row.names = FALSE)
write.table(data18, file = "/home/renshida/eqtm/data/human/30390659-cis.txt", sep = "\t", row.names = FALSE)
has_na <- any(is.na(data18$Genename))
na_count <- sum(is.na(data18$Genename))
print(has_na)
print(na_count)



#######
data19 <- read_excel("31791327-cis.xlsx")
data19 <- data19[-c(806:811),]
colnames(data19) <- data19[2,]
data19 <- data19[-c(1,2),]
data19.1 <- data19[,c(1,3,7,8,9)]
colnames(data19.1)[1] <- "CpG"
colnames(data19.1)[2] <- "Gene"
colnames(data19.1)[3] <- "Beta"
colnames(data19.1)[4] <- "SE"
colnames(data19.1)[5] <- "P-value"
data19.1$FDR <- NA
data19.1$Species <- "Human"
data19.1$Tissue <- "Whole blood"
data19.1$Disease <- "Never smokers"
data19.1$Article <- "DNA methylation is associated with lung function in never smokers"
data19.1$Journal <- "Respiratory Research"
data19.1$Date <- "2019.12"
data19.1$PMID <- "31791327"
data19.1$Cis_Trans_eqtm <- "Cis"
data19.1$Method <- "Robust linear regression"
data19.1$Sample <- "166"
a <- c(1,2,14,3:13,15,16)
data19.1 <- data19.1[,a]
write.table(data19.1, file = "/home/renshida/eqtm/data/human/31791327-cis-LL.txt", sep = "\t", row.names = FALSE)


data19.2 <- data19[,c(1,3,10,11,12)]
colnames(data19.2)[1] <- "CpG"
colnames(data19.2)[2] <- "Gene"
colnames(data19.2)[3] <- "Beta"
colnames(data19.2)[4] <- "SE"
colnames(data19.2)[5] <- "P-value"
data19.2$FDR <- NA
data19.2$Species <- "Human"
data19.2$Tissue <- "Whole blood"
data19.2$Disease <- "Never smokers"
data19.2$Article <- "DNA methylation is associated with lung function in never smokers"
data19.2$Journal <- "Respiratory Research"
data19.2$Date <- "2019.12"
data19.2$PMID <- "31791327"
data19.2$Cis_Trans_eqtm <- "Cis"
data19.2$Method <- "Robust linear regression"
data19.2$Sample <- "906"
a <- c(1,2,14,3:13,15,16)
data19.2 <- data19.2[,a]
write.table(data19.2, file = "/home/renshida/eqtm/data/human/31791327-cis-LLS.txt", sep = "\t", row.names = FALSE)

data19.3 <- data19[,c(1,3,13,14,15)]
colnames(data19.3)[1] <- "CpG"
colnames(data19.3)[2] <- "Gene"
colnames(data19.3)[3] <- "Beta"
colnames(data19.3)[4] <- "SE"
colnames(data19.3)[5] <- "P-value"
data19.3$FDR <- NA
data19.3$Species <- "Human"
data19.3$Tissue <- "Whole blood"
data19.3$Disease <- "Never smokers"
data19.3$Article <- "DNA methylation is associated with lung function in never smokers"
data19.3$Journal <- "Respiratory Research"
data19.3$Date <- "2019.12"
data19.3$PMID <- "31791327"
data19.3$Cis_Trans_eqtm <- "Cis"
data19.3$Method <- "Robust linear regression"
data19.3$Sample <- "206"
a <- c(1,2,14,3:13,15,16)
data19.3 <- data19.3[,a]
write.table(data19.3, file = "/home/renshida/eqtm/data/human/31791327-cis-NTR.txt", sep = "\t", row.names = FALSE)


data19.4 <- data19[,c(1,3,16,17,18)]
colnames(data19.4)[1] <- "CpG"
colnames(data19.4)[2] <- "Gene"
colnames(data19.4)[3] <- "Beta"
colnames(data19.4)[4] <- "SE"
colnames(data19.4)[5] <- "P-value"
data19.4$FDR <- NA
data19.4$Species <- "Human"
data19.4$Tissue <- "Whole blood"
data19.4$Disease <- "Never smokers"
data19.4$Article <- "DNA methylation is associated with lung function in never smokers"
data19.4$Journal <- "Respiratory Research"
data19.4$Date <- "2019.12"
data19.4$PMID <- "31791327"
data19.4$Cis_Trans_eqtm <- "Cis"
data19.4$Method <- "Robust linear regression"
data19.4$Sample <- "150"
a <- c(1,2,14,3:13,15,16)
data19.4 <- data19.4[,a]
write.table(data19.4, file = "/home/renshida/eqtm/data/human/31791327-cis-RS.txt", sep = "\t", row.names = FALSE)




#####
data20 <- read_tsv("27918535.txt")


data20 <- data20[,c(2,17,8,20,1,22)]
data20 <- data20 %>%
  separate(`Beta (SE)`, into = c("Beta1", "Beta2", "Beta3", "Beta4"), sep = ";\\s*", convert = TRUE)
data20 <- data20 %>%
  separate(Beta1, into = c("Beta1_value", "Beta1_SE"), sep = "\\(", convert = TRUE) %>%
  separate(Beta2, into = c("Beta2_value", "Beta2_SE"), sep = "\\(", convert = TRUE) %>%
  separate(Beta3, into = c("Beta3_value", "Beta3_SE"), sep = "\\(", convert = TRUE) %>%
  separate(Beta4, into = c("Beta4_value", "Beta4_SE"), sep = "\\(", convert = TRUE)
data20 <- data20 %>%
  mutate(across(everything(), ~ gsub("\\)", "", .)))
colnames(data20)[1] <- "CpG"
colnames(data20)[2] <- "Gene"
colnames(data20)[3] <- "Cis_Trans_eqtm"
data20RS <- data20[,c(1,2,3,4,5,12,13)]
colnames(data20RS)[4] <- "Beta"
colnames(data20RS)[5] <- "SE"
colnames(data20RS)[6] <- "P-value"
data20RS$Species <- "Human"
data20RS$Tissue <- "Whole blood"
data20RS$Disease <- "Chronic diseases"
data20RS$Article <- "Disease variants alter transcription factor levels and methylation of their binding sites"
data20RS$Journal <- "Nature Genetics"
data20RS$Date <- "2017.01"
data20RS$PMID <- "27918535"
data20RS$Method <- "Linear regression model"
data20RS$Sample <- "650"
data20RS$Cis_Trans_eqtm <- "Cis"
#range(data20RS$Beta)
#data20RS<- na.omit(data20RS)
#data20RS$Beta <- as.numeric(data20RS$Beta)
write.table(data20RS, file = "/home/renshida/eqtm/data/human/27918535-cis-RS.txt", sep = "\t", row.names = FALSE)

data20LLD <- data20[,c(1,2,3,6,7,12,13)]
colnames(data20LLD)[4] <- "Beta"
colnames(data20LLD)[5] <- "SE"
colnames(data20LLD)[6] <- "P-value"
data20LLD<- na.omit(data20LLD)
data20LLD$Species <- "Human"
data20LLD$Tissue <- "Whole blood"
data20LLD$Disease <- "Normal"
data20LLD$Article <- "Disease variants alter transcription factor levels and methylation of their binding sites"
data20LLD$Journal <- "Nature Genetics"
data20LLD$Date <- "2017.01"
data20LLD$PMID <- "27918535"
data20LLD$Method <- "Linear regression model"
data20LLD$Sample <- "625"
data20LLD$Cis_Trans_eqtm <- "Cis"
write.table(data20LLD, file = "/home/renshida/eqtm/data/human/27918535-cis-LLD.txt", sep = "\t", row.names = FALSE)


data20LLS <- data20[,c(1,2,3,8,9,12,13)]
colnames(data20LLD)[4] <- "Beta"
colnames(data20LLD)[5] <- "SE"
colnames(data20LLD)[6] <- "P-value"
data20LLS<- na.omit(data20LLS)
data20LLS$Species <- "Human"
data20LLS$Tissue <- "Whole blood"
data20LLS$Disease <- "Normal"
data20LLS$Article <- "Disease variants alter transcription factor levels and methylation of their binding sites"
data20LLS$Journal <- "Nature Genetics"
data20LLS$Date <- "2017.01"
data20LLS$PMID <- "27918535"
data20LLS$Method <- "Linear regression model"
data20LLS$Sample <- "652"
data20LLS$Cis_Trans_eqtm <- "Cis"
write.table(data20LLS, file = "/home/renshida/eqtm/data/human/27918535-cis-LLS.txt", sep = "\t", row.names = FALSE)


data20CO <- data20[,c(1,2,3,10,11,12,13)]
colnames(data20CO)[4] <- "Beta"
colnames(data20CO)[5] <- "SE"
colnames(data20CO)[6] <- "P-value"
data20CO<- na.omit(data20CO)
data20CO$Species <- "Human"
data20CO$Tissue <- "Whole blood"
data20CO$Disease <- "Type 2 diabetes/Cardiovascular disease"
data20CO$Article <- "Disease variants alter transcription factor levels and methylation of their binding sites"
data20CO$Journal <- "Nature Genetics"
data20CO$Date <- "2017.01"
data20CO$PMID <- "27918535"
data20CO$Method <- "Linear regression model"
data20CO$Sample <- "184"
data20CO$Cis_Trans_eqtm <- "Cis"
write.table(data20CO, file = "/home/renshida/eqtm/data/human/27918535-cis-CO.txt", sep = "\t", row.names = FALSE)






##########
data21 <- read_excel("35710981-cis.xlsx",sheet = "Table S12")
colnames(data21) <- data21[1,]
data21 <- data21[-1,]
data21 <- data21[,c(2,4,8,9,10)]
colnames(data21)[1] <- "CpG"
colnames(data21)[2] <- "Gene"
colnames(data21)[3] <- "Beta"
colnames(data21)[4] <- "P-value"
colnames(data21)[5] <- "FDR"
data21$Species <- "Human"
data21$Tissue <- "Kidney"
data21$Disease <- "Normal"
data21$Article <- "Epigenomic and transcriptomic analyses define core cell types, genes and targetable mechanisms for kidney disease"
data21$Journal <- "Nature Genetics"
data21$Date <- "2022.07"
data21$PMID <- "35710981"
data21$Method <- "Linear regression model"
data21$Sample <- 414
data21$Cis_Trans_eqtm <- "Cis"
data21$SE <- NA
a <- c(1:2,15,3:14,16)
data21 <- data21[,a]
a <- c(1:4,16,5:15)
data21 <- data21[,a]
write.table(data21, file = "/home/renshida/eqtm/data/human/35710981-cis.txt", sep = "\t", row.names = FALSE)



####
# data22 <- read_excel("33752734-cis-trans.xlsx",sheet = "Table S2")
# colnames(data22) <- data22[1,]
# data22 <- data22[-1,]
# data22 <- NULL



############
# data23 <- read_csv("33989148-cis.txt")
# str(data23)
# head(data23)
# data23 <- data23 %>%
  separate(`CpG Gene Statistic P-Value FDR Beta`, into = c("CpG", "Gene", "Statistic", "P.Value", "FDR", "Beta"), sep = " ",  convert = TRUE)
data23 <- data23[,c(1,2,6,4,5)]
colnames(data23)[1] <- "CpG"
colnames(data23)[2] <- "Gene"
colnames(data23)[3] <- "Beta"
colnames(data23)[4] <- "P-value"
colnames(data23)[5] <- "FDR"
data23$Species <- "Human"
data23$Tissue <- "Nasal epithelial"
data23$Disease <- "Differences between smokers and never-smokers among asthmatics"
data23$Article <- "Current Smoking Alters Gene Expression and DNA Methylation in the Nasal Epithelium of Patients with Asthma"
data23$Journal <- "American Journal of Respiratory Cell and Molecular Biology"
data23$Date <- "2021.10"
data23$PMID <- "33989148"
data23$Method <- "Linear regression model"
data23$Sample <- 29
data23$Cis_Trans_eqtm <- "Cis"
data23$SE <- NA
a <- c(1:3,16,4:15)
data23 <- data23[,a]
a <- c(1:2,16,3:15)
data23 <- data23[,a]
write.table(data23, file = "/home/renshida/eqtm/data/human/33989148-cis.txt", sep = "\t", row.names = FALSE)



##########
data24 <- read_excel("31578227-cis-DMR.xlsx",sheet = "Table S9")
colnames(data24) <- data24[1,]
data24 <- data24[-1,]
data24 <- data24[,c(1,6,17,13,15,16)]
colnames(data24)[1] <- "DMR"
colnames(data24)[2] <- "nCpG"
colnames(data24)[3] <- "Gene"
colnames(data24)[4] <- "Beta"
colnames(data24)[5] <- "P-value"
colnames(data24)[6] <- "FDR"
data24$Species <- "Human"
data24$Tissue <- "Pituitary"
data24$Disease <- "Growth Hormone–Secreting Pituitary Tumors"
data24$Article <- "Genetic and Epigenetic Characterization of Growth Hormone-Secreting Pituitary Tumors"
data24$Journal <- "Molecular Cancer Research"
data24$Date <- "2019.12"
data24$PMID <- "31578227"
data24$Method <- "Linear regression model"
data24$Sample <- 21
data24$Cis_Trans_eqtm <- "Cis"
data24$SE <- NA
a <- c(1:4,17,5:16)
data24 <- data24[,a]
a <- c(1:3,17,4:16)
data24 <- data24[,a]
write.table(data24, file = "/home/renshida/eqtm/data/human/31578227-cis-DMR.txt", sep = "\t", row.names = FALSE)

############
data25 <- read_excel("35232286-DMR.xlsx",sheet = "eQTM_P<0.05")
data25 <- data25[,c(5,10,1,6,11)]
colnames(data25)[1] <- "DMR"
colnames(data25)[2] <- "nCpG"
colnames(data25)[3] <- "Gene"
colnames(data25)[4] <- "Beta"
colnames(data25)[5] <- "P-value"
data25$FDR <- NA
data25$Species <- "Human"
data25$Tissue <- "Cerebrum"
data25$Disease <- "Down syndrome"
data25$Article <- "Prenatal NeuN+ neurons of Down syndrome display aberrant integrative DNA methylation and gene expression profiles"
data25$Journal <- "Epigenomics"
data25$Date <- "2022.04"
data25$PMID <- "35232286"
data25$Method <- "Pearson correlation coefficient"
data25$Sample <- 32
data25$Cis_Trans_eqtm <- NA
data25$SE <- NA
a <- c(1:4,17,5:16)
data25 <- data25[,a]
a <- c(1:3,17,4:16)
data25 <- data25[,a]
write.table(data25, file = "/home/renshida/eqtm/data/human/35232286-DMR.txt", sep = "\t", row.names = FALSE)


########
data26 <- read.table("36803404.txt", header = TRUE, sep = "\t")
#data26 <- separate(data26, col = 1, into = c("CpG", "Gene", "CpG_Chr", "Position_kb", "Transcript_Gene", "Coefficient", "p_value"), sep = "\t", convert = TRUE)
data26 <- data26[,c(2,5,6,7)]
colnames(data26)[1] <- "CpG"
colnames(data26)[2] <- "Gene"
colnames(data26)[3] <- "Beta"
colnames(data26)[4] <- "P-value"
data26$FDR <- NA
data26$Species <- "Human"
data26$Tissue <- "Peripheral blood"
data26$Disease <- "Systemic sclerosis"
data26$Article <- "Distinct genome-wide DNA methylation and gene expression signatures in classical monocytes from African American patients with systemic sclerosis"
data26$Journal <- "Clinical Epigenetics"
data26$Date <- "2023.02"
data26$PMID <- "36803404"
data26$Method <- "Linear regression model"
data26$Sample <- 34
data26$Cis_Trans_eqtm <- NA
data26$SE <- NA
a <- c(1:2,15,3:14,16)
data26 <- data26[,a]
a <- c(1:4,16,5:15)
data26 <- data26[,a]
write.table(data26, file = "/home/renshida/eqtm/data/human/36803404-top25.txt", sep = "\t", row.names = FALSE)




##########
data27 <- read_csv("32569636-cis-top30.txt")
head(data27)
data27 <- separate(data27, col = 1, into = c("Probe", "Gene", "Beta", "P", "FDR", "Beta_t", "P_t", "FDR_t"), sep = "\t", convert = TRUE)
data27e <- data27[,c(1,2,3,4,5)]
colnames(data27e)[1] <- "CpG"
colnames(data27e)[2] <- "Gene"
colnames(data27e)[3] <- "Beta"
colnames(data27e)[4] <- "P-value"
colnames(data27e)[5] <- "FDR"
data27e$Species <- "Human"
data27e$Tissue <- "Nasal epithelial"
data27e$Disease <- "Atopic asthma"
data27e$Article <- "Expression Quantitative Trait Methylation Analysis Reveals Methylomic Associations With Gene Expression in Childhood Asthma"
data27e$Journal <- "Chest"
data27e$Date <- "2020.11"
data27e$PMID <- "32569636"
data27e$Method <- "Multivariate linear regression"
data27e$Sample <- 455
data27e$Cis_Trans_eqtm <- "Cis"
data27e$SE <- NA
a <- c(1:3,16,4:15)
data27e <- data27e[,a]
a <- c(1:2,16,3:15)
data27e <- data27e[,a]
write.table(data27e, file = "/home/renshida/eqtm/data/human/32569636-cis-top30-eva.txt", sep = "\t", row.names = FALSE)

                                                                                                        
data27g <- data27[,c(1,2,6,7,8)]
colnames(data27g)[1] <- "CpG"
colnames(data27g)[2] <- "Gene"
colnames(data27g)[3] <- "Beta"
colnames(data27g)[4] <- "P-value"
colnames(data27g)[5] <- "FDR"
data27g$Species <- "Human"
data27g$Tissue <- "Nasal epithelial"
data27g$Disease <- "Atopic asthma"
data27g$Article <- "Expression Quantitative Trait Methylation Analysis Reveals Methylomic Associations With Gene Expression in Childhood Asthma"
data27g$Journal <- "Chest"
data27g$Date <- "2020.11"
data27g$PMID <- "32569636"
data27g$Method <- "Multivariate linear regression"
data27g$Sample <- 69
data27g$Cis_Trans_eqtm <- "Cis"
data27g$SE <- NA
a <- c(1:3,16,4:15)
data27g <- data27g[,a]
a <- c(1:2,16,3:15)
data27g <- data27g[,a]
write.table(data27g, file = "/home/renshida/eqtm/data/human/32569636-cis-top30-gse.txt", sep = "\t", row.names = FALSE)



########
data28 <- read_excel("37563237-cis-trans.xlsx",sheet = "ST1")
colnames(data28) <- data28[1,]
data28 <- data28[-1,]
data28 <- data28[,c(1,9,11,14)]
colnames(data28)[1] <- "CpG"
colnames(data28)[2] <- "Gene"
colnames(data28)[3] <- "Beta"
colnames(data28)[4] <- "P-value"
data28$FDR <- NA
data28$Species <- "Human"
data28$Tissue <- "Whole blood"
data28$Disease <- "Normal"
data28$Article <- "Expression quantitative trait methylation analysis elucidates gene regulatory effects of DNA methylation: the Framingham Heart Study"
data28$Journal <- "Scientific Reports"
data28$Date <- "2023.08"
data28$PMID <- "37563237"
data28$Method <- "Linear regression model"
data28$Sample <- 2115
data28$Cis_Trans_eqtm <- "Cis"
data28$SE <- NA
a <- c(1:3,16,4:15)
data28 <- data28[,a]
a <- c(1:2,16,3:15)
data28 <- data28[,a]
data28$`P-value` <- as.numeric(data28$`P-value`)
range(data28$`P-value`)
write.table(data28, file = "/home/renshida/eqtm/data/human/37563237-cis-top10000.txt", sep = "\t", row.names = FALSE)


data29 <- read_excel("37563237-cis-trans.xlsx",sheet = "ST2")
colnames(data29) <- data29[1,]
data29 <- data29[-1,]
data29 <- data29[,c(1,9,11,14)]
colnames(data29)[1] <- "CpG"
colnames(data29)[2] <- "Gene"
colnames(data29)[3] <- "Beta"
colnames(data29)[4] <- "P-value"
data29$FDR <- NA
data29$Species <- "Human"
data29$Tissue <- "Whole blood"
data29$Disease <- "Normal"
data29$Article <- "Expression quantitative trait methylation analysis elucidates gene regulatory effects of DNA methylation: the Framingham Heart Study"
data29$Journal <- "Scientific Reports"
data29$Date <- "2023.08"
data29$PMID <- "37563237"
data29$Method <- "Linear regression model"
data29$Sample <- 2115
data29$Cis_Trans_eqtm <- "Trans"
data29$SE <- NA
a <- c(1:3,16,4:15)
data29 <- data29[,a]
a <- c(1:2,16,3:15)
data29 <- data29[,a]
data29$`P-value` <- as.numeric(data29$`P-value`)
range(data29$`P-value`)
write.table(data29, file = "/home/renshida/eqtm/data/human/37563237-trans-top10000.txt", sep = "\t", row.names = FALSE)






