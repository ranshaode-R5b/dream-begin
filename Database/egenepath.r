BiocManager::install(c("ReactomePA"))
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(enrichplot)
install.packages("tidytext")
library(tidytext)
library(dplyr)
setwd("/home/renshida/eqtm/data/data_cpg/")
###31282290
data1 <- read.table("31282290-cis-mir.txt",header = TRUE,sep = "")
data1 <- data.frame(Gene = data1[, 2])
#data1$Gene <- str_replace(data1$Gene, "hsa-miR-(\\w+)(-\\d+p)?", "mir\\1")
gene_symbols <- unique(na.omit(data1$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
head(data1)
go_result <- enrichGO(gene          = gene_symbols,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "SYMBOL",
                      ont           = "BP",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)            # 将Entrez ID转换为基因符号
barplot(go_result, showCategory = 10, title = "GO Enrichment")

#####37271320
data2 <- read.table("37271320-cis.txt",header = TRUE,sep = "")
data2.1 <- read.table("37271320-trans1.txt",header = TRUE,sep = "")
data2.2 <- read.table("37271320-trans2.txt",header = TRUE,sep = "")
data2.3 <- read.table("37271320-trans3.txt",header = TRUE,sep = "")
data2 <- data.frame(Gene = data2[, 2])
data2.1 <- data.frame(Gene = data2.1[, 2])
data2.2 <- data.frame(Gene = data2.2[, 2])
data2.3 <- data.frame(Gene = data2.3[, 2])
data2.4<- rbind(data2,data2.1,data2.2,data2.3)
gene_symbols <- unique(na.omit(data2.4$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)            # 将Entrez ID转换为基因符号
go_df <- as.data.frame(go_result)
ontology_colors <- c("BP" = "darkcyan", "MF" = "sienna", "CC" = "steelblue")

# 绘制柱状分面图
ggplot(go_df %>% filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
         arrange(p.adjust) %>%
         group_by(ONTOLOGY) %>%
         slice_head(n = 10),
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(go_df %>% filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
         arrange(p.adjust) %>%
         group_by(ONTOLOGY) %>%
         slice_head(n = 10), aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +  # 确保在每个分面内正确排序
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()




##########39393618
data3 <- read.table("39393618-cis.txt",header = TRUE,sep="")
data3 <- data.frame(Gene = data3[, 2])
gene_symbols <- unique(na.omit(data3$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
go_df <- as.data.frame(go_result)
##柱状图
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


#######31076557
data4 <- read.table("31076557-cis.txt",header = TRUE,sep="")
data4 <- data.frame(Gene = data4[, 2])
gene_symbols <- unique(na.omit(data4$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


######35710981
data6 <-read.table("35710981-cis.txt",header = T,sep="")
data5 <- data.frame(Gene = data5[, 2])
gene_symbols <- unique(na.omit(data5$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


######31791327
data6 <-read.table("31791327-cis-LL.txt",header = T,sep="")
data6 <- data.frame(Gene = data6[, 2])
gene_symbols <- unique(na.omit(data6$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()




#####35302492
data7 <-read.table("35302492-cis-unadj.txt",header = T,sep="")
data7 <- data.frame(Gene = data7[, 2])
gene_symbols <- unique(na.omit(data7$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()



###########33092652
data8 <-read.table("33092652-cis-trans-kora.txt",header = T,sep="")
data8 <- data.frame(Gene = data8[, 2])
gene_symbols <- unique(na.omit(data8$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


####
data8.1 <-read.table("33092652-cis-trans-rot.txt",header = T,sep="")
data8.1 <- data.frame(Gene = data8.1[, 2])
gene_symbols <- unique(na.omit(data8.1$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


####29914364
data9 <-read.table("29914364-GTP.txt",header = T,sep="")
data9 <- data.frame(Gene = data9[, 2])
gene_symbols <- unique(na.omit(data9$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


###
data9.1 <-read.table("29914364-MESA.txt",header = T,sep="")
data9.1 <- data.frame(Gene = data9.1[, 2])
gene_symbols <- unique(na.omit(data9.1$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


######38438351
data10 <-read.table("38438351-cis.txt",header = T,sep="")
data10 <- data.frame(Gene = data10[, 2])
gene_symbols <- unique(na.omit(data10$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()



#####35426765
data11 <-read.table("38438351-cis.txt",header = T,sep="")
data11 <- data.frame(Gene = data11[, 2])
gene_symbols <- unique(na.omit(data11$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


#####37233989
data12 <-read.table("37233989-cis-vc.txt",header = T,sep="")
data12.1 <- read.table("37233989-cis-ve.txt",header = T,sep="")
data12 <- data.frame(Gene = data12[, 2])
data12.1 <- data.frame(Gene = data12.1[, 2])
data12 <- rbind(data12,data12.1)
gene_symbols <- unique(na.omit(data12$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()


######37598461
data13 <-read.table("37598461-cis.txt",header = T,sep="")
data13 <- data.frame(Gene = data13[, 2])
gene_symbols <- unique(na.omit(data13$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()



getwd()
####32928872
data14 <-read.table("32928872-cis-cite.txt",header = T,sep="")
data14 <- data.frame(Gene = data14[, 2])
gene_symbols <- unique(na.omit(data14$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()



####32493484
data15 <-read.table("32493484-cis.txt",header = T,sep="")
data15 <- data.frame(Gene = data15[, 2])
gene_symbols <- unique(na.omit(data15$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()





#######31578227
data16 <-read.table("31578227-cis-DMR.txt",header = T,sep="")
data16 <- data.frame(Gene = data16[, 3])
gene_symbols <- unique(na.omit(data16$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()




#####35232286
data17 <-read.table("35232286-DMR.txt",header = T,sep="")
data17 <- data.frame(Gene = data17[, 3])
gene_symbols <- unique(na.omit(data17$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()




#######30390659
data18 <-read.table("30390659-cis.txt",header = T,sep="")
data18 <- data.frame(Gene = data18[, 2])
gene_symbols <- unique(na.omit(data18$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()



#####37563237
data19 <-read.table("37563237-cis-top10000.txt",header = T,sep="")
data19.1 <-read.table("37563237-trans-top10000.txt",header = T,sep="")
data19 <- data.frame(Gene = data19[, 2])
data19.1 <- data.frame(Gene = data19[, 2])
data19 <- rbind(data19,data19.1)
gene_symbols <- unique(na.omit(data19$Gene))
gene.df <- bitr(gene_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
any(is.na(gene.df))
go_result <- enrichGO(gene          = gene.df$ENTREZID,
                      OrgDb         = org.Hs.eg.db,
                      keyType       = "ENTREZID",
                      ont           = "ALL",            # 选择 "BP", "MF" 或 "CC"
                      pAdjustMethod = "BH",
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)    

go_df <- as.data.frame(go_result)
##柱状图
top_terms <- go_df %>%
  filter(ONTOLOGY %in% c("BP", "MF", "CC")) %>%
  arrange(p.adjust) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>%
  ungroup() %>%
  mutate(Description = reorder_within(Description, -p.adjust, ONTOLOGY))
ggplot(top_terms,
       aes(x = reorder(Description, -p.adjust), y = -log10(p.adjust), fill = ONTOLOGY)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  facet_wrap(~ ONTOLOGY, scales = "free", ncol = 1) +
  scale_x_reordered() +
  scale_fill_manual(values = ontology_colors) +
  theme_minimal() +
  labs(
    x = "GO Term",
    y = "-log10(p.adjust)",
    fill = "Ontology",
    title = "Top 10 GO Terms Enrichment"
  ) +
  theme(axis.text = element_text(face = "bold", color = "gray50"))

####气泡图
ggplot(top_terms, aes(x = -log10(p.adjust), y = Description, size = Count, fill = p.adjust)) +
  geom_point(shape = 21, color = 'black') +
  facet_grid(ONTOLOGY ~ ., scales = 'free_y', space = 'free_y') +
  scale_fill_gradient(low = '#E27371', high = '#5D82A7') +
  scale_y_reordered() +
  labs(title = 'GO Enrichment',
       y = 'GO term',
       x = '-log10(Adjusted p-value)') +
  guides(fill = guide_colorbar(reverse = TRUE, order = 1)) +
  theme_bw()





