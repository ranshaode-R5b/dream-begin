#安装seurat包
#BiocManager::install("Seurat")
#BiocManager::install("harmony")
#install.packages(c('dplyr','patchwork','Rcpp','ggplot2'))

#加载包
library(Seurat)
library(Rcpp)
library(harmony)
library(dplyr)
library(patchwork)
library(ggplot2)
setwd("home/renshida/scrna-Sle")
getwd()
dir()

#数据导入
#control
data_dir <- "./control/"
list.files(data_dir) 
con <- Read10X(data.dir = data_dir)
dim(con)   #即基因数和细胞数
con[1:10,1:4] 
#查看数据矩阵中线粒体基因表达情况
con[grep(rownames(con),pattern="^MT-"),1:3]
rowSums(con[grep(rownames(con),pattern="^MT-"),])


#treat
data_dir <- "./stim/"
list.files(data_dir)
stim <- Read10X(data.dir = data_dir) 
dim(stim)


save(con,stim,file = "sample.Rda")
load("sample.Rda") 
#样本合并
con <- CreateSeuratObject(counts = con, project = "control", min.cells = 3, min.features = 200)
stim <- CreateSeuratObject(counts = stim, project = "stimulus", min.cells = 3, min.features = 200)
seurat_object <- merge(con,stim,add.cell.ids = c("con","stim"))


#计算线粒体基因百分比，存放到 metadata 中； 
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
VlnPlot(seurat_object, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#ggsave("plot1.pdf", width = 28, height = 25, units = "cm")
#筛选数据 
seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 5)

##绘制散点图，查看相关性
plot1 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat_object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2


#数据标准化：LogNormalize 的算法：A = log( 1 + ( UMIA ÷ UMITotal ) × 10000 ) 
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000)
#鉴定细胞间表达量高变的基因，用于下游分析，PCA
seurat_object <- FindVariableFeatures(seurat_object, selection.method = "vst", nfeatures = 2000)
# 提取表达量变变化最高的 10 个基因； 
top10 <- head(VariableFeatures(seurat_object), 10)
top10
# 绘制带有和不带有标签的变量特征的散点图
plot3 <- VariableFeaturePlot(seurat_object)+NoLegend()
plot4 <- LabelPoints(plot = plot3, points = top10, repel = TRUE, xnudge=0, ynudge=0)
plot3+plot4

#排除细胞周期基因对分群和降维的影响
cc.genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
#计算细胞周期分数
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])#在Seurat v5中，GetAssayData函数无法处理包含多个层（layers）的assay
seurat_object <- CellCycleScoring(seurat_object, s.features = s.genes, g2m.features = g2m.genes, set.ident = FALSE)
head(seurat_object[[]])
#去除细胞周期对分群和降维的影响
seurat_object <- ScaleData(seurat_object, vars.to.regress = c("S.Score", "G2M.Score"), features = VariableFeatures(seurat_object))
#而对所有基因进行标准化的方法如下： 
#seurat_object <- ScaleData(seurat_object, vars.to.regress = c("S.Score", "G2M.Score"), features = row.names(seurat_object)) ##耗时2min
#线性降维（PCA）,默认用高变基因集
seurat_object <- RunPCA(seurat_object, features = VariableFeatures(object = seurat_object)) 
#检查 PCA 分群结果， 这里只展示前 12 个 PC,每个 PC 只显示 3 个基因； 
print(seurat_object[["pca"]], dims = 1:12, nfeatures = 3) 
DimPlot(seurat_object, reduction = "pca")
#画前 2 个主成分的热图
DimHeatmap(seurat_object, dims = 1:2, cells = 500, balanced = TRUE) 
#ggsave("plot2.png", width = 15, height = 13, units = "cm")



#Jackstraw 置换检验算法；重复取样（原数据的 1%），重跑PCA,鉴定p-value较小的PC；计算‘null distribution’(即零假设成立时)时的基因 scores; 
seurat_object <- JackStraw(seurat_object, num.replicate = 100,dims = 20)  ##耗时3min
seurat_object <- ScoreJackStraw(seurat_object, dims = 1:20) #dims必须小于等于上一步
JackStrawPlot(seurat_object, dims = 1:15) #dims必须在上一步的范围内
#肘部图，基于每个主成分对方差解释率的排名； 
ElbowPlot(seurat_object)
#使用拐点作为降维的度数
seurat_object <- FindNeighbors(seurat_object, dims = 1:15)
seurat_object <- FindClusters(seurat_object, resolution = 0.5) 
#使用 Idents（）函数可查看不同细胞的分群； 查看前5个细胞的分群ID
head(Idents(seurat_object), 5) 

#tsne非线性降维
seurat_object <- RunTSNE(seurat_object, dims = 1:15) 
#用TSNEPlot函数绘制tsne图
tsneplot<-TSNEPlot(seurat_object,label = TRUE, pt.size = 1.5)+ NoLegend()
tsneplot
#用DimPlot函数绘制分组tsne图
tsneplot1 <- DimPlot(seurat_object, reduction = "tsne", group.by = "orig.ident", pt.size = 1.5)
tsneplot1
#绘制 Marker 基因的 tsne 图； 
FeaturePlot(seurat_object, features = c("MS4A1", "CD14")) 
#UMAP非线性降维
seurat_object <- RunUMAP(seurat_object, dims = 1:15)
head(seurat_object@reductions$umap@cell.embeddings)
#用DimPlot函数绘制UMAP图
umapplot <- DimPlot(seurat_object, reduction = "umap",pt.size = 1.5,label = T)
umapplot


#为分群重新指定细胞类型 
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8+ T","FCGR3A+ Mono", "NK", "DC", "Platelet","T","Eryth","Mk","HSPC","CD8- T")
names(new.cluster.ids) 
levels(seurat_object)
#将seurat_object的水平属性赋值给new.cluster.ids的names属性； 
names(new.cluster.ids) <- levels(seurat_object)
names(new.cluster.ids) 
seurat_object <- RenameIdents(seurat_object, new.cluster.ids)  
tsneplot2<-TSNEPlot(seurat_object,label = TRUE, pt.size = 1.4)
tsneplot2
#ggsave("tsne.png",tsneplot2,width = 20, height = 15, units = "cm")
#保存工作空间
save(seurat_object,file = "uncorrected_obj.Rda")



#亚群间差异基因分析
dif<-FindAllMarkers(seurat_object,logfc.threshold = 0.6,min.pct = 0.4,only.pos = T)
#logfc.threshold定义上调倍数阈值，min.pct定义基因至少在细胞亚群中多少细胞中表达，only.pos确定只筛选上调基因
sig.dif<-dif%>%group_by(cluster)%>%top_n(n = 5,wt = avg_log2FC)
#保存差异分析结果
write.table(sig.dif,"sig.dif.xls",row.names = T,col.names = NA,quote = F,sep = "\t")

#组间差异基因分析
#首先为细胞添加新的标签，以同时区分分组和细胞亚群
seurat_object[["group.cluster"]]<-paste(seurat_object$orig.ident,Idents(seurat_object),sep = '_')
#比较B细胞在两组中的差异基因
group.dif<-FindMarkers(seurat_object,group.by = "group.cluster",
                       ident.1 = "control_B",ident.2 = "stimulus_B",
                       logfc.threshold = 0.6,min.pct = 0.4)
#差异倍数最大的top10差异基因
group.sig.dif<-group.dif%>%top_n(n = 10,wt = abs(avg_log2FC))
#上调倍数最大和最小的各10个基因
group.sig.dif1<-rbind(group.dif%>%top_n(n = 10,wt = avg_log2FC),group.dif%>%top_n(n = -10,wt = avg_log2FC))
write.table(group.sig.dif,"group.sig.dif.xls",row.names = T,col.names = NA,quote = F,sep = "\t")





#批次效应矫正#
###Harmony矫正
rm(list = ls())  #清除内存所有变量
load("sample.Rda")  #加载数据
#数据集中测到的少于200个基因的细胞（min.features = 200）和少于3个细胞覆盖的基因（min.cells = 3）被过滤掉
con <- CreateSeuratObject(counts = con, project = "control", min.cells = 3, min.features = 200)
stim <- CreateSeuratObject(counts = stim, project = "stimulus", min.cells = 3, min.features = 200)
#合并为一个seurat对象，主要是为了添加样本标识和过滤细胞
seurat_object <- merge(con,stim,add.cell.ids = c("con","stim"))
seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object,pattern = "^MT-")
seurat_object <- subset(seurat_object,nFeature_RNA > 200 & nFeature_RNA < 1500 & percent.mt < 5)
#归一化、高变基因筛选、标准化和降维
seurat_object <- NormalizeData(seurat_object, normalization.method = "LogNormalize", scale.factor = 10000) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData()%>%RunPCA()
#harmony矫正
seurat_object = RunHarmony(seurat_object,"orig.ident", plot_convergence = TRUE)
#细胞聚类
seurat_object <- seurat_object %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.5)
#降维可视化
p1 <- DimPlot(seurat_object, reduction = "umap", group.by = "orig.ident", pt.size = 0.5)
p2 <- DimPlot(seurat_object, reduction = "umap", pt.size = 0.5)
p1+p2

#为分群重新指定细胞类型 
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8+ T","FCGR3A+ Mono", "NK", "DC", "Platelet","T","Eryth","Mk","HSPC","CD8- T") #自定义名称
names(new.cluster.ids)
levels(seurat_object)
names(new.cluster.ids) <- levels(seurat_object)
names(new.cluster.ids)
seurat_object <- RenameIdents(seurat_object, new.cluster.ids) 
P1<-UMAPPlot(seurat_object,label = TRUE, pt.size = 1.4) 
P1
save(seurat_object,file = "harmony_obj.Rda")
