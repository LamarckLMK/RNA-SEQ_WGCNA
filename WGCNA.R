library(WGCNA)
library(DESeq2)

# enableWGCNAThreads(nThreads = 10)

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
dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-28],colData = sample_metadata,design = ~ Zone)  # 创建DESeqDataSet对象
mcols(dataExpr_deseq)$basepairs = data0$geneLengt1  # 为DESeqDataSet添加基因长度信息
fpkm_matrix = fpkm(dataExpr_deseq)  # 计算FPKM矩阵
datExpr = t(log2(fpkm_matrix+1))  # 对FPKM矩阵进行对数转换

# head(datExpr[1:5,1:5]) # samples in row, genes in column
# match(sample_metadata$sample_ID, colnames(data0))
# datExpr <- datExpr[,1:1000]

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

# 画一条水平线  选择软阈值时  R^最好要≥0.85
abline(h=0.85,col="red")

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


################################################################################
################################################################################
# Option 2a: step-by-step
power = power
adjacency = adjacency(datExpr, power = power)
TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
dissTOM = 1-TOM

#===============================================================================
#
#  Construct modules (proceed with the genetree from option 2b)
#
#===============================================================================
# Plot gene tree
geneTree = hclust(as.dist(dissTOM), method = "average");
#pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);
dev.off()

# Module identification using dynamic tree cut
# We like large modules, so we set the minimum module size relatively high:
# minModuleSize = 30;
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
                            pamRespectsDendro = FALSE,minClusterSize = 30);
table(dynamicMods)
length(table(dynamicMods)) 
# Convert numeric labels into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
pdf(file = "C:/Users/Lamarck/Desktop/module_tree.pdf", width = 8, height = 6);
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
                    hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
dev.off()

#===============================================================================
#
#  Merge modules
#
#===============================================================================
# Calculate eigengenes
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")

# Merge close modules
MEDissThres=0.40
abline(h=MEDissThres, col = "red")
merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
mergedColors = merge$colors  
mergedMEs = merge$newMEs  
# Plot merged module tree
#pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
                    c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
dev.off()



#===============================================================================
#
#  PART 1: Correlate module eigen-genes and samples (or other discrete data)
#
#===============================================================================

# 加载必要的包
library("pheatmap")

# 第一个热图：旧模块特征基因热图
pdf(file = "C:/Users/Lamarck/Desktop/oldMEs_heatmap.pdf", width = 8, height = 6)
rownames(merge$oldMEs) = names(data0[, -28])  # 设置行名
pheatmap(merge$oldMEs,
         cluster_col = TRUE, cluster_row = TRUE,
         show_rownames = FALSE, show_colnames = TRUE,
         fontsize = 6)
dev.off()  # 结束PDF保存

# 准备数据和注释信息
col_ann <- sample_metadata[, c(1, 3)]  # 提取样本注释信息
rownames(col_ann) <- col_ann[, 1]  # 设置注释信息行名
col_ann <- data.frame(col_ann)
col_ann$Zone <- as.factor(col_ann$Zone)  # 将分组信息转换为因子
col_ann <- col_ann[order(col_ann$Zone), ]  # 按Zone排序
col_ann$sample_ID <- NULL  # 删除多余列
ann_color <- list(
  "Zone" = c("Z1" = "blue", "Z2" = "green")  # 设置分组颜色
)

data <- data.frame(merge$newMEs)  # 提取新模块特征基因数据
data <- data[order(match(rownames(data), rownames(col_ann))), ]  # 排序
rownames(merge$newMEs) = names(data0[, -28])  # 设置行名

# 第二个热图：新模块特征基因与样本特征的热图
pdf(file = "C:/Users/Lamarck/Desktop/newMEs_heatmap.pdf", width = 8, height = 6)
pheatmap(data,
         cluster_col = TRUE, cluster_row = FALSE,
         show_rownames = FALSE, show_colnames = TRUE,
         fontsize = 6,
         annotation_row = col_ann, annotation_colors = ann_color)
dev.off()  # 结束PDF保存

