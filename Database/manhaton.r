setwd("/home/renshida/eqtm/data/data_cpg/")
install.packages("qqman")
library(qqman)
install.packages("CMplot")
library(CMplot)
########31282290
data1 <- read.table("31282290-cis-mir.txt",header = TRUE,sep = "")
data1 <- data.frame(data1[, c(1,19,20,6)])
colnames(data1)[1:4] <- c("SNP","CHR","BP","P")
data1$CHR <- gsub("chr", "", data1$CHR)
data1$CHR <- gsub("X", "23", data1$CHR)
data1$CHR <- gsub("Y", "24", data1$CHR)
data1$CHR <- as.numeric(data1$CHR)
#manhattan(data1, main = "eQTM manhattan", col = c("blue4", "orange3"), genomewideline = 3, ylim = c(0, 10))
table(data1$CHR)

# 绘制曼哈顿图
#CMplot(data1, plot.type = "m", threshold = 0.01, threshold.col = c('grey', 'black'),chr.border = TRUE,
#       threshold.lty = c(1, 2), threshold.lwd = c(1, 1), amplify = TRUE,
#       signal.cex = c(1, 1), signal.pch = c(20, 20), signal.col = colors,
#       col = colors, file.output = FALSE)
colors <- c("blue", "green", "purple")
CMplot(data1, plot.type = "m", col = colors, threshold = 0.01, threshold.col = 'black',chr.border = TRUE,
       threshold.lty = 2, threshold.lwd = 1,amplify = TRUE,height = 8, width = 20, dpi = 300)


#####37271320
# data2 <- read.table("37271320-cis.txt",header = TRUE,sep = "")
# data2.1 <- read.table("37271320-trans1.txt",header = TRUE,sep = "")
# data2.2 <- read.table("37271320-trans2.txt",header = TRUE,sep = "")
# data2.3 <- read.table("37271320-trans3.txt",header = TRUE,sep = "")
# data2.4<- rbind(data2,data2.1,data2.2,data2.3)
# write.table(data2.4, file = "/home/renshida/eqtm/data/data_cpg/37271320-all.txt", sep = "\t", row.names = FALSE)
data2 <- read.table("37271320-all.txt", header = TRUE, sep = "")
data2 <- data.frame(data2[, c(1,19,20,7)])
colnames(data2)[1:4] <- c("SNP","CHR","BP","P")
data2$CHR <- gsub("chr", "", data2$CHR)
data2$CHR <- gsub("X", "23", data2$CHR)
data2$CHR <- gsub("Y", "24", data2$CHR)
data2$CHR <- as.numeric(data2$CHR)
#manhattan(data2, main = "eQTM manhattan", col = c("blue4", "orange3"), genomewideline = -log10(5e-2), ylim = c(0, 30))
range(data2$P)
table(data2$CHR)
CMplot(data2,plot.type = "m",threshold = 0.05, threshold.col='black',
                 threshold.lty = 1,threshold.lwd = 1, amplify = T,
                 signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c("red","orange"))
CMplot(data2, plot.type = "m", col = colors, threshold = 0.05, threshold.col = 'black',chr.border = TRUE,
                 threshold.lty = 2, threshold.lwd = 1,amplify = TRUE,height = 8, width = 20, dpi = 300)
str(data2$CHR)
table(data2$CHR)
######39393618
data3 <- read.table("39393618-cis.txt",header = TRUE,sep="")
data3 <- data.frame(data3[, c(1,19,20,7)])
colnames(data3)[1:4] <- c("SNP","CHR","BP","P")
data3$CHR <- gsub("chr", "", data3$CHR)
data3$CHR <- gsub("X", "23", data3$CHR)
data3$CHR <- gsub("Y", "24", data3$CHR)
data3$CHR <- as.numeric(data3$CHR)
#manhattan(data2, main = "eQTM manhattan", col = c("blue4", "orange3"), genomewideline = -log10(5e-2), ylim = c(0, 30))
range(data3$P)

CMplot(data3,plot.type = "m",threshold = 0.05, threshold.col='black',
       threshold.lty = 1,threshold.lwd = 1, amplify = T,
       signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c("red","orange"),file.output = FALSE)

CMplot(data3, plot.type = "m", col = colors, threshold = 0.05, threshold.col = 'black',chr.border = TRUE,
       threshold.lty = 2, threshold.lwd = 1,amplify = TRUE,height = 8, width = 20, dpi = 300)

table(data3$CHR)
###31791327
data4 <- read.table("31791327-cis-LL.txt",header=TRUE,sep="")
data4 <- data.frame(data4[, c(1,19,20,6)])
range(data4$P.value)
colnames(data4)[1:4] <- c("SNP","CHR","BP","P")
data4$CHR <- gsub("chr", "", data4$CHR)
data4$CHR <- gsub("X", "23", data4$CHR)
data4$CHR <- gsub("Y", "24", data4$CHR)
data4$CHR <- as.numeric(data4$CHR)
#manhattan(data2, main = "eQTM manhattan", col = c("blue4", "orange3"), genomewideline = -log10(5e-2), ylim = c(0, 30))
range(data4$P)

CMplot(data4,plot.type = "m",threshold = 0.05, threshold.col='black', chr.border = TRUE,
       threshold.lty = 1,threshold.lwd = 1, amplify = T,
       signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c("red","orange"),height = 8, width = 20, dpi = 300)
table(data4$CHR)




###31076557
data5 <- read.table("31076557-cis.txt",header=T,sep="")
data5 <- data.frame(data5[, c(1,19,20,7)])
colnames(data5)[1:4] <- c("SNP","CHR","BP","P")
data5$CHR <- gsub("chr", "", data5$CHR)
data5$CHR <- gsub("X", "23", data5$CHR)
data5$CHR <- gsub("Y", "24", data5$CHR)
data5$CHR <- as.numeric(data5$CHR)
#manhattan(data2, main = "eQTM manhattan", col = c("blue4", "orange3"), genomewideline = -log10(5e-2), ylim = c(0, 30))
range(data5$P)
data5$CHR <- as.numeric(as.character(data5$CHR))
#CMplot(data5,plot.type = "m",threshold = 0.01, threshold.col='black', chr.border = TRUE,
       #threshold.lty = 1,threshold.lwd = 1, amplify = T,
       #signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c("red","orange"))
CMplot(data5,plot.type = "m", col = colors, threshold = 0.01, threshold.col = 'black',chr.border = TRUE,
       threshold.lty = 2, threshold.lwd = 1,amplify = TRUE,height = 8, width = 20, dpi = 300)
table(data5$CHR)
head(data5)
# data5$CHR <- as.factor(data5$CHR)
# # 计算全基因组上的累积坐标
# chr_offsets <- data5 %>%
#   group_by(CHR) %>%
#   summarise(chr_min = min(BP)) %>%
#   mutate(offset = cumsum(lag(chr_min, default = 0)))  # 计算染色体偏移量
# 
# # 合并累积坐标
# data5 <- data5 %>%
#   left_join(chr_offsets, by = "CHR") %>%
#   mutate(cum_BP = BP + offset)
# 
# # 获取染色体中点位置
# chr_midpoints <- data5 %>%
#   group_by(CHR) %>%
#   summarise(mid = mean(cum_BP))
# 
# # 修改横坐标，使用 `scale_x_discrete` 进行更好的控制
# ggplot(data5, aes(x = as.factor(CHR), y = -log10(P), color = as.factor(CHR))) +
#   geom_point(alpha = 0.7, size = 1) +  # 调整点的透明度和大小
#   geom_hline(yintercept = -log10(0.01), linetype = "dashed", color = "black", linewidth = 0.5) +  # 添加显著性阈值
#   scale_x_discrete(
#     breaks = as.character(chr_midpoints$CHR),  # 确保染色体标签有合理间距
#     labels = as.character(chr_midpoints$CHR)  # 标签使用染色体号
#   ) +
#   scale_color_manual(values = color_palette) +  # 使用颜色
#   labs(x = "Chromosome", y = "-log10(P-value)", title = "Manhattan Plot") +
#   theme_minimal() +
#   theme(
#     legend.position = "none",  # 不显示图例
#     panel.grid.major.x = element_line(color = "grey80", linetype = "dashed"),  # 在染色体之间加虚线
#     axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
# )


#####35710981
data6 <- read.table("35710981-cis.txt",header=T,sep="")
data6<- data.frame(data6[, c(1,19,20,7)])
colnames(data6)[1:4] <- c("SNP","CHR","BP","P")
data6$CHR <- gsub("chr", "", data6$CHR)
data6$CHR <- gsub("X", "23", data6$CHR)
data6$CHR <- gsub("Y", "24", data6$CHR)
data6$CHR <- as.numeric(data6$CHR)
#manhattan(data2, main = "eQTM manhattan", col = c("blue4", "orange3"), genomewideline = -log10(5e-2), ylim = c(0, 30))
range(data6$P)
data6$CHR <- as.numeric(as.character(data6$CHR))
#CMplot(data5,plot.type = "m",threshold = 0.01, threshold.col='black', chr.border = TRUE,
#threshold.lty = 1,threshold.lwd = 1, amplify = T,
#signal.cex = c(1,1), signal.pch = c(20,20),signal.col = c("red","orange"))
CMplot(data6,plot.type = "m", col = colors, threshold = 0.05, threshold.col = 'black',chr.border = TRUE,
       threshold.lty = 2, threshold.lwd = 1,amplify = TRUE,height = 8, width = 20, dpi = 300)
table(data5$CHR)



















































































