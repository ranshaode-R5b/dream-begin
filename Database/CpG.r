sessionInfo()

##########33092652
#BiocManager::install("illuminaHumanv4.db")

# 加载所需的包
library(illuminaHumanv4.db)
library(AnnotationDbi)

# 假设 combined_data 数据框中，Gene 列存储了 Illumina 探针 ID
probe_ids <- combined_data$Gene

# 利用 mapIds 函数将 Illumina 探针 ID 转换为 Gene Symbol
gene_symbols <- mapIds(illuminaHumanv4.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# 将转换结果添加为新的列
combined_data$GeneSymbol <- gene_symbols
has_na <- any(is.na(combined_data$GeneSymbol))
na_count <- sum(is.na(combined_data$GeneSymbol))
print(has_na)
print(na_count)
write.table(combined_data[,c(2,17)], file = "/home/renshida/eqtm/data/trans_g_id/333092652-cis-trans-rot.txt", sep = "\t", row.names = FALSE)


#####combined_data1
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
write.table(combined_data1[,c(2,17)], file = "/home/renshida/eqtm/data/trans_g_id/33092652-cis-trans-kora.txt", sep = "\t", row.names = FALSE)




######38438351
gene_ids <- data13$Gene
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = gene_ids,
                       column = "SYMBOL",
                       keytype = "ENSEMBL",
                       multiVals = "first")
# 将基因符号添加到 data13 数据框
data13$GeneSymbol <- gene_symbols
has_na <- any(is.na(data13$GeneSymbol))
na_count <- sum(is.na(data13$GeneSymbol))
print(has_na)
print(na_count)
write.table(data13[,c(2,17)], file = "/home/renshida/eqtm/data/trans_g_id/38438351-cis.txt", sep = "\t", row.names = FALSE)



#####31076557
library(AnnotationDbi)
library(org.Hs.eg.db)
data15$Gene <- sub("\\..*", "", data15$Gene)
data15$GeneSymbol <- mapIds(org.Hs.eg.db,
                            keys = data15$Gene,
                            column = "SYMBOL",
                            keytype = "ENSEMBL",
                            multiVals = "first")
has_na <- any(is.na(data15$GeneSymbol))
na_count <- sum(is.na(data15$GeneSymbol))
print(has_na)
print(na_count)
write.table(data15[,c(2,17)], file = "/home/renshida/eqtm/data/trans_g_id/31076557-cis.txt", sep = "\t", row.names = FALSE)




#######cpg id 
# BiocManager::install("AnnotationDbi")
# BiocManager::install("org.Hs.eg.db")
# library(AnnotationDbi)
# library(org.Hs.eg.db)
# 假设你有CpG的位置信息


library(readr)
library(stringr)
# 读取数据，并将第七行作为列名
data <- read.csv("humanmethylation450_15017482_v1-2.csv", skip = 7, stringsAsFactors = FALSE)


############32928872   hg19版本cpg id匹配
data14_parsed <- data14 %>%
  mutate(
    Chromosome = str_extract(CpG, "chr[0-9XY]+"),
    Position = as.numeric(str_extract(CpG, "(?<=:)[0-9]+"))
  )
data$CHR <- as.character(data$CHR)
data$MAPINFO <- as.character(data$MAPINFO)

# 检查 data 中的连接列是否唯一
any(duplicated(data[c("CHR", "MAPINFO")]))
any(is.na(data14_parsed[c("Chromosome", "Position")]))

# 检查 data 中 CHR 和 MAPINFO 的缺失值
any(is.na(data[c("CHR", "MAPINFO")]))
data <- data[!duplicated(data[c("CHR", "MAPINFO")]), ]

# 移除缺失值
data14_parsed <- data14_parsed[!is.na(data14_parsed$Chromosome) & !is.na(data14_parsed$Position), ]
data <- data[!is.na(data$CHR) & !is.na(data$MAPINFO), ]

# 确保数据类型一致性
data14_parsed$Position <- as.character(data14_parsed$Position)
data$MAPINFO <- as.character(data$MAPINFO)

# 数据匹配
data14_parsed$Chromosome <- gsub("^chr", "", data14_parsed$Chromosome)
matched_data <- data14_parsed %>%
  left_join(data, by = c("Chromosome" = "CHR", "Position" = "MAPINFO"))

# 将匹配到的探针 ID 添加到 data14 中
data14_parsed$CpG_ID <- matched_data$IlmnID  # 假设探针 ID 列名为 "IlmnID"
has_na <- any(is.na(data14_parsed$CpG_ID))
na_count <- sum(is.na(data14_parsed$CpG_ID))
print(has_na)
print(na_count)





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

are_columns_equal <- all(ann$IlmnID == ann$Name)
# 输出结果
if (are_columns_equal) {
  print("两列信息完全一致")
} else {
  print("两列信息不一致")
}







#######hg38版本的cpg id匹配
# 确保 data 中的 CHR 和 MAPINFO 列存在，并转换为字符类型
data <- ann
data$hg38_chr <- as.character(data$hg38_chr)
data$hg38_pos <- as.character(data$hg38_pos)
# 确保数据类型一致性
data14_parsed$Position <- as.character(data14_parsed$Position)
data14_parsed$Chromosome <- as.character(data14_parsed$Chromosome)
# 数据匹配
data14_parsed$Chromosome <- gsub("^chr", "", data14_parsed$Chromosome)
data$hg38_chr <- gsub("^chr", "", data$hg38_chr)
matched_data <- data14_parsed %>%
  left_join(data, by = c("Chromosome" = "hg38_chr", "Position" = "hg38_pos"))

# 将匹配到的探针 ID 添加到 data14 中
data14_parsed$CpG_ID <- matched_data$IlmnID  # 假设探针 ID 列名为 "IlmnID"




#########31403346 cpgid转换
library(stringr)

split_result <- str_split(data1$CpG, "[.][.]+")

# 提取染色体和位置信息
data1$Chromosome <- sapply(split_result, function(x) {
  # 提取染色体编号，去除 "chr" 前缀
  chrom <- strsplit(x[1], "\\.")[[1]][1]
  gsub("chr", "", chrom)
})
data1$Position <- sapply(split_result, function(x) {
  # 提取位置信息，取第一个数字部分
  pos <- strsplit(x[2], "\\.")[[1]][1]
  pos
})
#hg19
data$MAPINFO <- as.character(data$MAPINFO)
matched_data <- data1 %>%
  left_join(data, by = c("Chromosome" = "CHR", "Position" = "MAPINFO"))

# 将匹配到的探针 ID 添加到 data14 中
data1$CpG_ID <- matched_data$IlmnID  # 假设探针 ID 列名为 "IlmnID"
#hg38
data <- ann
data$hg38_pos <- as.character(data$hg38_pos)
matched_data <- data1 %>%
  left_join(data, by = c("Chromosome" = "hg38_chr", "Position" = "hg38_pos"))

# 将匹配到的探针 ID 添加到 data14 中
data1$CpG_ID <- matched_data$IlmnID  # 假设探针 ID 列名为 "IlmnID"






######根据gtfv38提取gene symbol  31076557
gtf <- read.csv("gene_id_name.txt",header = FALSE, sep = "")
colnames(gtf)[1] <- "ens"
colnames(gtf)[2] <- "sym"
head(gtf)
data15$Gene <- sub("\\..*", "", data15$Gene)
gtf$ens <- sub("\\..*", "", gtf$ens)
data15 <- data15 %>%
  left_join(gtf, by = c("Gene" = "ens"))
# data15 <- data15[,c(-17,-18)]
has_na <- any(is.na(data15$sym))
na_count <- sum(is.na(data15$sym))
print(has_na)
print(na_count)




#####38438351
data13$Gene <- sub("\\..*", "", data13$Gene)
data13 <- data13 %>%
  left_join(gtf, by = c("Gene" = "ens"))
# data15 <- data15[,c(-17,-18)]
has_na <- any(is.na(data13$sym))
na_count <- sum(is.na(data13$sym))
print(has_na)
print(na_count)

########
#####

