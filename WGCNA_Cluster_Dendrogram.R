library(WGCNA)
library(DESeq2)

# 在处理数据框(data.frame)时，不会自动给将String类型转换成factor类型
options(stringsAsFactors = FALSE);

#===============================================================================
#
#  STEP1: Read the gene counts table and plot the sample tree
#
#===============================================================================

# 读取基因的表达矩阵
data0=read.table("C:/Users/Lamarck/Desktop/mycounts_lengths.txt",header=T,row.names=1,sep="\t")

# 读取分组文件
sample_metadata = read.csv(file = "C:/Users/Lamarck/Desktop/sample_info.csv")

# 计算表达矩阵中的fpkm值
dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-28],colData = sample_metadata,design = ~ class)  # 创建DESeqDataSet对象
mcols(dataExpr_deseq)$basepairs = data0$geneLengt1  # 为DESeqDataSet添加基因长度信息
fpkm_matrix = fpkm(dataExpr_deseq)  # 计算FPKM矩阵
datExpr = t(log2(fpkm_matrix+1))  # 对FPKM矩阵进行对数转换

# 计算sample distance  得到sample tree  看看有没有outlier
sampleTree = hclust(dist(datExpr), method = "average");

# 绘制sample tree
pdf(file = "C:/Users/Lamarck/Desktop/sample_tree.pdf", width = 40, height = 9);  #设置PDF页面的宽度和高度
par(cex = 1.3);  # 图形中所有文本元素大小的放大倍率
par(mar = c(0,4,2,0))  # 设置图形四个边的边距
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",
     cex.lab = 1.5,cex.axis = 1.5, cex.main = 2)  # 绘制层次聚类树
dev.off()

#===============================================================================
#
#  Step2: Choose soft threshold parameter
#
#===============================================================================

# 选择一系列的软阈值参数
powers = c(c(1:20), seq(from = 22, to=30, by=2))

# 计算不同软阈值下的网络拓扑拟合指数
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 设置图形参数  1行2列  同一张画布上绘制两个图
par(mfrow = c(1,2))

# 设置文本标签的大小
cex1 = 0.9;

# 绘制第一个图  Scale independence  标度独立性
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));  # 绘制标度独立性图
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");  # 将文本标签添加到图中

# 画一条水平线  选择软阈值时  R^最好要≥0.90
abline(h=0.90,col="red")

# 画第二个图  Mean connectivity  平均连接性  表示软阈值和平均连接性的关系
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))  # 绘制平均连接性图
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")  # 把文本标签添加到图中

dev.off()  # 结束绘图

#===============================================================================
#
#  STEP3: Turn data expression into topological overlap matrix
#
#===============================================================================

# 自动评估一个合适的软阈值
power=sft$powerEstimate

# Option 1: automatic  自动将基因表达数据转化为拓朴重叠矩阵TOM  并使用该矩阵进行模块识别
cor <- WGCNA::cor  # 设置相关系数函数

# 进行模块识别
net = blockwiseModules(datExpr, power = power,
                       TOMType = "signed", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = FALSE,
                       verbose = 3)
cor<- stats::cor  #恢复cor函数

# 绘制模块树状图  
sizeGrWindow(12, 9)  # 设置绘图窗口的大小
mergedColors = labels2colors(net$colors)  # 将模块标签转换为颜色
pdf(file = "C:/Users/Lamarck/Desktop/Cluster_Dendrogram.pdf", width = 8, height = 6);  # 保存绘图到PDF文件
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)  #  绘制模块树状图和模块颜色
dev.off()
