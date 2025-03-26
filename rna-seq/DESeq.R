setwd("~/rna-seq/data/subread")
library(DESeq2)
library(ggplot2)
#install.packages("pheatmap")
#remotes::install_version("enrichplot")
#install.packages("https://cran.r-project.org/src/contrib/Archive/scatterpie/scatterpie_0.1.8.tar.gz", repos=NULL)
#install.packages("clusterProfiler")
#install.packages("RColorBrewer")
#BiocManager::install("org.Rn.eg.db")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("topGO") 
library(RColorBrewer)
library(pheatmap)
library(clusterProfiler)
library(dplyr)
library(org.Rn.eg.db)
#library(org.Hs.eg.db)
getwd()
countdata <- read.table("all.id.txt", header = TRUE, skip = 1, row.names = 1)
countdata <- countdata[ ,c(-1:-5)]
head(countdata)
sampleNames <- c("SRR27175483", "SRR27175484", "SRR27175485", "SRR27175486", "SRR27175487", "SRR27175488")
names(countdata)[1:6] <- sampleNames
metadata <- read.csv("SraRunTable.csv", row.names = 1)
metadata$sampleid <- row.names(metadata)
metadata <- metadata[match(colnames(countdata), metadata$sampleid), ]
metadata <- metadata[,32:33]
group <- rep(c("TMAO", "control"), each=3)
ddsMat <- DESeqDataSetFromMatrix(countData = countdata,
                                 colData = metadata,
                                 design = ~treatment)
ddsMat <- DESeq(ddsMat)
results <- results(ddsMat, pAdjustMethod = "fdr", alpha = 0.05)
summary(results)
resultsNames(ddsMat)


# 添加基因全名
results$description <- mapIds(x = org.Rn.eg.db,
                              keys = row.names(results),
                              column = "GENENAME",
                              keytype = "SYMBOL",
                              multiVals = "first")

# 添加基因 symbol
results$symbol <- row.names(results)

# 添加 ENTREZ ID
results$entrez <- mapIds(x = org.Rn.eg.db,
                         keys = row.names(results),
                         column = "ENTREZID",
                         keytype = "SYMBOL",
                         multiVals = "first")

# 添加 ENSEMBL
results$ensembl <- mapIds(x = org.Rn.eg.db,
                          keys = row.names(results),
                          column = "ENSEMBL",
                          keytype = "SYMBOL",
                          multiVals = "first")
results <- as.data.frame(results)
# 取 (q < 0.05) 的基因
results_sig <- subset(results, padj < 0.05)

# 查看结果
head(results_sig)  


#热图
# 将所有样本转换为 rlog
ddsMat_rlog <- rlog(ddsMat, blind = FALSE)

# 收集40个显著基因，制作矩阵
mat <- assay(ddsMat_rlog[row.names(results_sig)])[1:40, ]

# 选择要用来注释列的列变量
annotation_col <- data.frame(
  Group = factor(colData(ddsMat_rlog)$treatment), 
  row.names = rownames(colData(ddsMat_rlog))
)
z <- factor(colData(ddsMat_rlog)$treatment) 
# 指定要用来注释列的颜色
ann_colors = list(
  Group = c(TMAO = "lightblue", control = "darkorange")
)

# 使用 pheatmap 功能制作热图
pheatmap(mat = mat, 
         color = colorRampPalette(brewer.pal(9, "YlOrBr"))(255), 
         scale = "row",
         annotation_col = annotation_col, 
         annotation_colors = ann_colors,
         fontsize = 6.5, 
         cellwidth = 55,
         show_colnames = F)



#火山图
# 从 DESeq2 结果中收集倍数变化和 FDR 校正的 pvalue
## - 将 pvalues 更改为 -log10 (1.3 = 0.05)
data <- data.frame(gene = row.names(results),
                   pval = -log10(results$padj), 
                   lfc = results$log2FoldChange)

# 删除任何以 NA 的行
data <- na.omit(data)
data <- mutate(data, color = case_when(data$lfc > 2 & data$pval > 1.3 ~ "Increased",
                                       data$lfc < 2 & data$pval > 1.3 ~ "Decreased",
                                       data$pval < 1.3 ~ "nonsignificant"))
table(data$color)
# 用 x-y 值制作一个基本的 ggplot2 对象
vol <- ggplot(data, aes(x = lfc, y = pval, color = color))

# 添加 ggplot2 图层
vol +   
  ggtitle(label = "Volcano Plot", subtitle = "Colored by fold-change direction") +
  geom_point(size = 2.5, alpha = 0.8, na.rm = T) +
  scale_color_manual(name = "Directionality",
                     values = c(Increased = "#008B00", Decreased = "#CD4F39", nonsignificant = "darkgray")) +
  theme_bw(base_size = 14) + 
  theme(legend.position = "right") + 
  xlab(expression(log[2]("TMAO" / "control"))) + 
  ylab(expression(-log[10]("adjusted p-value"))) + 
  geom_hline(yintercept = 1.3, colour = "darkgrey") + 
  scale_y_continuous(trans = "log1p") 



#MA
library(ggplot2)
# 将padj值转换为-log10以便于绘图
results$neg_log10_padj <- -log10(results$padj)
results$logbasemean <- log10(results$baseMean)

results_sig <- as.data.frame(results_sig)
BiocManager::install("ggrepel")
library(ggrepel)
results <- as.data.frame(results)
significant_genes <- results %>% 
  filter(padj < 0.05) %>% 
  head(20)
ggplot(results, aes(x = logbasemean, y = log2FoldChange)) +
  geom_point(aes(color = neg_log10_padj), alpha = 0.6, size = 2) +  
  scale_color_gradient(low = "blue", high = "red") +  # 颜色渐变从蓝色到红色
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # 添加y=0的虚线
  geom_vline(xintercept = c(0), linetype = "dashed", color = "grey") +  # 添加x=0的虚线
  geom_text_repel(data = significant_genes, aes(label = symbol), size = 3, color = "black") +  # 标记显著基因
  labs(x = "Average Expression Level (logbasemean)", 
       y = "Log2 Fold Change (M)", 
       color = "Padj") +
  theme_minimal()  # 使用简洁主题
#2
ggplot(results, aes(x = logbasemean, y = log2FoldChange)) +
  geom_point(aes(color = "Significant"), alpha = 0.6, size = 2) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  # 添加y=0的虚线
  geom_vline(xintercept = c(0), linetype = "dashed", color = "grey") +  # 添加x=0的虚线
  geom_point(data = significant_genes, aes(color = "Significant"), alpha = 0.6, size = 2) +  # 单独标记显著基因
  geom_text_repel(data = significant_genes, aes(label = symbol), size = 3, color = "black") +  # 标记显著基因
  scale_color_manual(values = c("Significant" = "red")) +  # 为显著基因指定颜色
  labs(x = "Average Expression Level (logbasemean)", 
       y = "Log2 Fold Change (M)", 
       color = "Significance") +
  theme_minimal()



results <- as.data.frame(results)

# 为上调和下调的基因指定颜色
results$color <- ifelse(results$log2FoldChange > 0 & data$neg_log10_padj > 1.3, "Upregulated", "Downregulated")
results$color <- ifelse(results$log2FoldChange > 0 & results$neg_log10_padj > 1.3, "Significant Upregulated", 
                        ifelse(results$log2FoldChange < 0 & results$neg_log10_padj > 1.3, "Significant Downregulated", "Not Significant"))

# 筛选出padj值小于0.05的基因作为显著基因
significant_genes <- results %>% 
  filter((results$log2FoldChange > 0 & results$neg_log10_padj > 1.3) | 
           (results$log2FoldChange < 0 & results$neg_log10_padj > 1.3)) %>% 
  head(20)

# 绘制MA图，使用不同颜色表示上调和下调的基因
ggplot(results, aes(x = logbasemean, y = log2FoldChange, color = color)) +
  geom_point(alpha = 0.6, size = 2) +  
  scale_color_manual(values = c("Significant Upregulated" = "red", "Significant Downregulated" = "blue", "Not Significant" = "grey")) +  
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +  
  geom_vline(xintercept = c(0), linetype = "dashed", color = "grey") +  
  geom_text_repel(data = significant_genes, aes(label = symbol), size = 3, color = "black") +  
  labs(x = "Average Expression Level (logbasemean)", 
       y = "Log2 Fold Change (M)", 
       color = "Regulation") +
  theme_minimal()
#基因通路富集


# 删除没有任何 entrez 标识符的基因
results_sig_entrez <- subset(results_sig, is.na(entrez) == FALSE)
# 创建一个log2倍数变化的基因矩阵
gene_matrix <- results_sig_entrez$log2FoldChange
# 添加 entrezID 作为每个 logFC 条目的名称
names(gene_matrix) <- results_sig_entrez$entrez
# 查看基因矩阵的格式
##- Names = ENTREZ ID
##- Values = Log2 Fold changes
head(gene_matrix)  



#GO
go_enrich <- enrichGO(gene = names(gene_matrix),
                      OrgDb = 'org.Rn.eg.db', 
                      readable = T,
                      ont = "BP",
                      pvalueCutoff = 0.05, 
                      qvalueCutoff = 0.10)
# 结果可视化
barplot(go_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "GO Biological Pathways",
        font.size = 8)



#KEGG
kegg_enrich <- enrichKEGG(gene = results_sig_entrez$entrez,
                          keyType = "kegg",
                          organism = 'rno',
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10)

# 结果可视化
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)





df_results <- read.csv("df_results.csv",row.names = 1)
new_df_results <- df_results[!sapply(df_results$entrez,is.na),]
kegg_enrich <- enrichKEGG(gene = new_df_results$entrez,
                          keyType = "kegg",
                          organism = "mmu",
                          pvalueCutoff = 0.05, 
                          qvalueCutoff = 0.10,
                          )

# 结果可视化
barplot(kegg_enrich, 
        drop = TRUE, 
        showCategory = 10, 
        title = "KEGG Enrichment Pathways",
        font.size = 8)
