library(WGCNA)
library(DESeq2)

# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);  # 在处理数据框(data.frame)时，不会自动给将String类型转换成factor类型

#===============================================================================
#
#  STEP1: Read the gene counts table and plot the sample tree
#
#===============================================================================

# 读取基因的表达矩阵
data0=read.table("C:/Users/Lamarck/Desktop/gene_counts_table_WGCNA_LC.txt",header=T,row.names=1,sep="\t")

# 读取分组文件
sample_metadata = read.csv(file = "C:/Users/Lamarck/Desktop/sample_info.csv")

# 计算表达矩阵中的fpkm值
dataExpr_deseq <- DESeqDataSetFromMatrix(countData = data0[,-181],colData = sample_metadata,design = ~ Zone)  # 创建DESeqDataSet对象
mcols(dataExpr_deseq)$basepairs = data0$geneLengt1  # 为DESeqDataSet添加基因长度信息
fpkm_matrix = fpkm(dataExpr_deseq)  # 计算FPKM矩阵
datExpr = t(log2(fpkm_matrix+1))  # 对FPKM矩阵进行对数转换

# 查看fpkm矩阵以及截取前5000个基因
head(datExpr[1:5,1:5])  # 查看datExpr矩阵的前5行前5列
match(sample_metadata$sample_ID, colnames(data0))  # 查找sample_metadata中sample_ID在data0列名中的位置
datExpr <- datExpr[,1:5000]  # 截取前5000个基因的数据

# 计算sample distance  得到sample tree  看看有没有outlier
sampleTree = hclust(dist(datExpr), method = "average");

# 绘制sample tree
pdf(file = "1-n-sampleClustering.pdf", width = 40, height = 9);  #设置PDF页面的宽度和高度
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

# unsigned -> nodes with positive & negative correlation are treated equally 
# signed -> nodes with negative correlation are considered *unconnected*, treated as zero

# 绘制模块树状图  
sizeGrWindow(12, 9)  # 设置绘图窗口的大小
mergedColors = labels2colors(net$colors)  # 将模块标签转换为颜色
pdf(file = "4-module_tree_blockwise.pdf", width = 8, height = 6);  # 保存绘图到PDF文件
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)  #  绘制模块树状图和模块颜色
dev.off()

########################################################################################

# # Option 2a: step-by-step
# power = power
# adjacency = adjacency(datExpr, power = power)
# TOM = TOMsimilarity(adjacency); # Turn adjacency into topological overlap
# dissTOM = 1-TOM
# 
# 
# # Option 2b: 
# TOM = TOMsimilarityFromExpr(datExpr, power = power)
# dissTOM = 1-TOM 
# dim(dissTOM)
# 
# #===============================================================================
# #
# #  Construct modules (proceed with the genetree from option 2b)
# #
# #===============================================================================
# 
# # Plot gene tree
# geneTree = hclust(as.dist(dissTOM), method = "average");
# #pdf(file = "3-gene_cluster.pdf", width = 12, height = 9);
# plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
#      labels = FALSE, hang = 0.04);
# dev.off()
# 
# # Module identification using dynamic tree cut
# # We like large modules, so we set the minimum module size relatively high:
# # minModuleSize = 30;
# dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,deepSplit = 2, 
#                             pamRespectsDendro = FALSE,minClusterSize = 30);
# table(dynamicMods)
# length(table(dynamicMods)) 
# # Convert numeric labels into colors
# dynamicColors = labels2colors(dynamicMods)
# table(dynamicColors)
# # Plot the dendrogram and colors underneath
# pdf(file = "4-module_tree.pdf", width = 8, height = 6);
# plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",dendroLabels = FALSE,
#                     hang = 0.03,addGuide = TRUE, guideHang = 0.05,main = "Gene dendrogram and module colors")
# dev.off()
# 
# #===============================================================================
# #
# #  Merge modules
# #
# #===============================================================================
# 
# # Calculate eigengenes
# MEList = moduleEigengenes(datExpr, colors = dynamicColors)
# MEs = MEList$eigengenes
# # Calculate dissimilarity of module eigengenes
# MEDiss = 1-cor(MEs);
# # Cluster module eigengenes
# METree = hclust(as.dist(MEDiss), method = "average");
# # Plot the result
# sizeGrWindow(7, 6)
# plot(METree, main = "Clustering of module eigengenes",
#      xlab = "", sub = "")
# 
# # Merge close modules
# MEDissThres=0.40
# abline(h=MEDissThres, col = "red")
# merge = mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3) 
# mergedColors = merge$colors  
# mergedMEs = merge$newMEs  
# # Plot merged module tree
# #pdf(file = "5-merged_Module_Tree.pdf", width = 12, height = 9)  
# plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), 
#                     c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, 
#                     hang = 0.03, addGuide = TRUE, guideHang = 0.05)  
# dev.off()
# #write.table(merge$oldMEs,file="oldMEs.txt");
# #write.table(merge$newMEs,file="newMEs.txt");

########################################################################################


#===============================================================================
#
#  Export of networks to external software
#
#===============================================================================

# Export the gene list of old modules 
for (i in 1:length(merge$oldMEs)){
  modules = c(substring(names(merge$oldMEs)[i], 3));
  genes = colnames(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/orign_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/orign_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}
# Export the gene list of new modules 
for (i in 1:length(merge$newMEs)){
  modules = c(substring(names(merge$newMEs)[i], 3));
  genes = colnames(datExpr)
  inModule = is.finite(match(dynamicColors,modules))
  modGenes = genes[inModule]
  modTOM=TOM[inModule,inModule]
  dimnames(modTOM)=list(modGenes,modGenes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("output_for_cytoscape/merge_CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                                 nodeFile = paste("output_for_cytoscape/merge_CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                                 weighted = TRUE, threshold = -1, nodeNames = modGenes, nodeAttr = dynamicColors[inModule]);
}

#===============================================================================
#
#  PART 1: Correlate module eigen-genes and samples (or other discrete data)
#
#===============================================================================
# Heatmap of old module eigen-genes and samples
#pdf(file="oldMEs.pdf",heigh=80,width=20)
library("pheatmap")
rownames(merge$oldMEs)=names(data0[,-181])
pheatmap(merge$oldMEs,cluster_col=T,cluster_row=T,show_rownames=F,show_colnames=T,fontsize=6)
dev.off()


# Heatmap of new module eigen-genes and sample trait (e.g. Zone)
col_ann <- sample_metadata[,c(1,3)]
rownames(col_ann) <- col_ann[,1]
col_ann <- data.frame(col_ann)
col_ann$Zone <- as.factor(col_ann$Zone)
col_ann <- col_ann[order(col_ann$Zone),]
col_ann$sample_ID <- NULL
head(col_ann)
ann_color <- list("col_ann" = c("Z1" = "yellow",
                                "Z2" = "red",
                                "Z3" = "green"))

data <- data.frame(merge$newMEs)
data <- data[order(match(rownames(data), rownames(col_ann))),]
dim(merge$newMEs)

#pdf(file="newMEs.pdf",heigh=60,width=20)
rownames(merge$newMEs)=names(data0[,-181])
pheatmap(data,cluster_col=T,cluster_row=F,show_rownames=F,
         show_colnames=T,fontsize=6,
         annotation_row = col_ann, annotation_colors = ann_color)
dev.off()

#=====================================================================================
#
#  PART 2: Correlation between gene modules and microbial traits (continuous data)
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, mergedColors)$eigengenes
MEs = orderMEs(MEs0)

# Read microbial data as traits
bac_traits = read.table("b_order_234.txt", header = T, sep = "\t")
rownames(bac_traits) = bac_traits[, 1]
bac_traits = bac_traits[, -1]
# sample names should be consistent in eigen genes and traits !!!!
bac_traits = bac_traits[match(rownames(MEs), rownames(bac_traits)), ]
table(rownames(MEs) == rownames(bac_traits))

# Calculate pearson correlation coefficients between module eigen-genes and traits
moduleTraitCor = cor(MEs, bac_traits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#write.table(moduleTraitCor,file="moduleTrait_correlation.txt");
#write.table(moduleTraitPvalue,file="moduleTrait_pValue.txt");


#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
pdf("5-module-traits-bacteria-order.pdf", width = 80, height = 15)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(bac_traits),
               yLabels = colnames(MEs),
               ySymbols = colnames(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#  Plot heatmap of module-traits relationship
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor[1:25,1:25], 2), "\n(",
                    signif(moduleTraitPvalue[1:25,1:25], 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor[1:25,1:25])
pdf("5-module-traits-bacteria-order1.pdf", width = 20, height = 10)
par(mar = c(15, 12, 5, 5));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor[1:25,1:25],
               xLabels = colnames(bac_traits[1:25,1:25]),
               yLabels = colnames(MEs[1:25,1:25]),
               ySymbols = colnames(MEs[1:25,1:25]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()

#=====================================================================================
#
#   Intramodular analysis: identifying genes with high geneModuleMembership & geneTraitSignficance
#
#=====================================================================================


# Define variable Verru containing the Verrucomicrobiales column of bac_traits
Verru = as.data.frame(bac_traits$Verrucomicrobiales);
names(Verru) = "Verrucomicrobiales"

# names (colors) of the modules
modNames = substring(names(MEs), 3)

MET = orderMEs(cbind(MEs, Verru))

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, Verru, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Verru), sep="");
names(GSPvalue) = paste("p.GS.", names(Verru), sep="");

# Plot the dendrogram
sizeGrWindow(6,6);
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(10,10,1,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


module = "lightgreen"
# Rename to moduleColors
moduleColors = mergedColors
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for Verrucomicrobiales",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "blue")


## Draw bubble plot for particular module
colsum_bac_traits <- colSums(bac_traits)
colsum_bac_traits <- data.frame(colsum_bac_traits)
colsum_bac_traits$b_order <- rownames(colsum_bac_traits)
library(tidyr)
moduleTraitCor_long <- data.frame(moduleTraitCor)
moduleTraitCor_long$module <- rownames(moduleTraitCor)
moduleTraitCor_long <- moduleTraitCor_long[,c(235,1:234)]
moduleTraitCor_long <- gather(moduleTraitCor_long, b_order, PCC, Pseudomonadales:Others, factor_key = TRUE)

moduleTraitPvalue_long <- data.frame(moduleTraitPvalue)
moduleTraitPvalue_long$module <- rownames(moduleTraitPvalue)
moduleTraitPvalue_long <- moduleTraitPvalue_long[,c(235,1:234)]
moduleTraitPvalue_long <- gather(moduleTraitPvalue_long, b_order, pval, Pseudomonadales:Others, factor_key = TRUE)

moduleTrait_long <- merge(moduleTraitCor_long, moduleTraitPvalue_long, by = c("module","b_order"))

bubble_Data <- merge(moduleTrait_long, colsum_bac_traits, by = "b_order")
#just want module = "lightgreen"
bubble_Data_lightgreen <- bubble_Data[which(bubble_Data$module == "MElightgreen"),]

library(ggplot2)
ggplot(bubble_Data_lightgreen, aes(x= colsum_bac_traits, y= PCC, size = colsum_bac_traits,
                                   color = PCC, label = b_order)) +
  geom_text(hjust = 1, size=3) +
  geom_point(alpha=1) + ylab("Module-taxon correlation") + xlab("Relative abundance (sum)") +
  theme_bw()


############# Summary ###################################

head(datExpr)[1:5,1:5] # transcriptome data

head(sample_metadata)[1:5,] # metadata (sample info)
head(bac_traits)[1:5,1:5] # external trait






#=====================================================================================
#
#   Cytoscape
#
#=====================================================================================


#if(!"RCy3" %in% installed.packages()){
#  install.packages("BiocManager")
#  BiocManager::install("RCy3")
#}

# https://cytoscape.org/cytoscape-automation/for-scripters/R/notebooks/
library(RCy3)

cytoscapePing () # make sure cytoscape is open
cytoscapeVersionInfo ()

###### for yellow module of the merged data (newMEs) #################################
edge <- read.delim("output_for_cytoscape/merge_CytoscapeInput-edges-lightgreen.txt")
colnames(edge)
colnames(edge) <- c("source", "target","weight","direction","fromAltName","toAltName")

node <- read.delim("output_for_cytoscape/merge_CytoscapeInput-nodes-lightgreen.txt")
colnames(node)  
colnames(node) <- c("id","altName","node_attributes") 

createNetworkFromDataFrames(node,edge[1:50,], title="my first network", collection="DataFrame Example")

################ customise the network visualization ##################################
# use other pre-set visual style
setVisualStyle('Marquee')

# set up my own style
style.name = "myStyle"
defaults <- list(NODE_SHAPE="diamond",
                 NODE_SIZE=30,
                 EDGE_TRANSPARENCY=120,
                 NODE_LABEL_POSITION="W,E,c,0.00,0.00")
nodeLabels <- mapVisualProperty('node label','id','p')
nodeFills <- mapVisualProperty('node fill color','node_attributes','d',c("A","B"), c("#FF9900","#66AAAA"))
arrowShapes <- mapVisualProperty('Edge Target Arrow Shape','interaction','d',c("activates","inhibits","interacts"),c("Arrow","T","None"))
edgeWidth <- mapVisualProperty('edge width','weight','p')

createVisualStyle(style.name, defaults, list(nodeLabels,nodeFills,arrowShapes,edgeWidth))
setVisualStyle(style.name)
