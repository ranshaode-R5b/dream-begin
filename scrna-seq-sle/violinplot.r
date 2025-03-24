#加载相关R包
library(Seurat)
library(ggplot2)

#读取数据
load('harmony_obj.Rda')
dim(seurat_object)

#从seurat对象提取基因的表达量信息
#提取基因IL4I1的表达量信息
seurat_object[["celltype"]]<-Idents(seurat_object)
IL4I1<-FetchData(seurat_object,vars = 'IL4I1')

#提取T细胞marker基因平均表达量信息
Tcell.marker<-c('CD3D','CD3E','CD3G')
Tcell.data<-seurat_object@assays$RNA@data[Tcell.marker,]
Tcell.data<-expm1(Tcell.data)
T_cell<-log1p(colMeans(Tcell.data))

#构建绘制小提琴图所需的表达量表格
volin.data<-data.frame(cell.id = colnames(seurat_object),
                       celltype = seurat_object$celltype,
                       IL4I1,T_cell)

#使用数据绘制小提琴图
color<-c("#00468B","#925E9F","#759EDD","#0099B4","#76D1B1","#42B540","#B8D24D",
         "#EDE447","#FAB158","#FF7777","#FD0000","#AD002A","#AE8691","#CE9573",
         "#756455")
p1<-ggplot(data = volin.data,aes(x = celltype,y = T_cell))+
  geom_violin(aes(fill = celltype),alpha = 0.8,scale = 'count')+
  scale_fill_manual(values = color)
p1
#在小提琴图的图形上增加离散点
p2<-p1+geom_jitter(alpha = 0.3,col = "black",show.legend = FALSE,width = 0.1,size = 2,pch = 20)
p2

#修改坐标轴名称和图片名称
p3<-p2+labs(x="Cell Type",y="Expression Levels",title = "T cell genes")
p3

#自定义主题
mytheme<-theme(panel.background = element_rect(fill = 'white',color = 'black'),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_line(color = 'lightgrey',size = 0.8),
               axis.text.y = element_text(size = 12,color = 'black'),
               axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1,size = 12,color = 'black'),
               axis.title = element_text(face = 'bold',color = 'black',size = 14),
               plot.title = element_text(size = 16,face = 'bold',color = 'black',hjust = 0.5),
               legend.title = element_blank(),
               legend.text = element_text(size = 14,color = 'black'))
p4<-p3+mytheme
p4
#将小提琴图图层放置在散点图之上，凸显小提琴图
p5<-ggplot(volin.data,aes(x = celltype,y = T_cell))+
  geom_jitter(col = 'black',show.legend = FALSE,width = 0.1,size = 1,pch = 20)+
  geom_violin(aes(fill = celltype),alpha = 0.6,scale = 'count')+
  scale_fill_manual(values = color)+
  labs(x="Cell Type",y="Expression Levels",title = "T cell genes")+
  mytheme
p5

#保存图片
ggsave(filename = "violin.pdf",plot = p5)
