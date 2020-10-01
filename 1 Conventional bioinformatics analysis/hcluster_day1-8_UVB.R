# 执行前设置====================================
# 清空暂存数据
rm(list=ls())
# 载入R包
library(pheatmap)
# 设置工作目录
setwd("E:\\联合利华\\L-DNB_3Dmodel_data\\L-DNB for UVB_SBT ctrl-ref\\PCA visualization based on DNB\\heatmap")

# 整理数据集====================================
# 载入数据
dataset1 <- read.table('UVB_all_sample_data.csv',header = TRUE, row.names = 1, sep = ',')
dataset2 <- read.table('UVB_no_day8-uvb_data.csv',header = TRUE, row.names = 1, sep = ',')
# 截取表达矩阵的一部分数据来绘制热图
exp_ds = dataset1
#  exp_ds = dataset[c(1:60),]
# 构建样本分类数据
Group=c(rep('UVB',25))
Time=c(rep(1,7),
       rep(2,2),
       rep(3,2),
       rep(4,3),
       rep(5,2),
       rep(6,7),
       rep(7,2),
       rep(8,3)

       )

annotation_c <- data.frame(Group, Time)
rownames(annotation_c) <- colnames(exp_ds)

# 绘制热图=====================================
pheatmap(exp_ds, #表达数据
         cluster_rows = T,#行聚类
         cluster_cols = T,#列聚类
         annotation_col =annotation_c, #样本分类数据
         annotation_legend=TRUE, # 显示样本分类
         show_rownames = F,# 显示行名
         show_colnames = T,# 显示列名
         scale = "row", #对行标准化
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100), # 热图基准颜色
         angle_col = 45, # 设置列偏转角度，可选 270, 0, 45, 90, 315，
         gaps_row = T,
         clustering_method = "average",
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation"
         
)

datExpr = as.data.frame(t(exp_ds))
# hclusts聚类算法, dist计算基因之间的距离
sampleTree = hclust(dist(datExpr), method = "average")
plot(sampleTree, main = "Sample clustering", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)





