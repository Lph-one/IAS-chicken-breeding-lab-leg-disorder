## 在CA和CT组中，分别进行WGCNA
setwd("/Users/liupeihao/Desktop")

install.packages("BiocManager")
BiocManager::install("WGCNA")
library(WGCNA)
options(stringsAsFactors = FALSE)

## CA组的WGCNA
# 读入代谢物数据
metaData <- read.csv("CAmeta.csv", stringsAsFactors = FALSE)
datExpr0 <- as.data.frame(t(metaData[, -1]))
names(datExpr0) <- metaData$ID
rownames(datExpr0) <- names(metaData)[-1]

# 检查样本与代谢物是否合格
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) {
  if(sum(!gsg$goodGenes) > 0)
    print(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ",")))
  if(sum(!gsg$goodSamples) > 0)
    print(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ",")))
  datExpr0 <- datExpr0[gsg$goodSamples, gsg$goodGenes]}
## 再次确认
gsg <- goodSamplesGenes(datExpr0, verbose = 3)
if (!gsg$allOK) stop("Data still not clean")

# 样本聚类，检测离群样本
sampleTree <- hclust(dist(datExpr0), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", xlab = "", sub = "")

# 按高度切割去除异常样本（可调整cutHeight）
clust <- cutreeStatic(sampleTree, cutHeight = 8e6, minSize = 10)
table(clust)
keepSamples <- (clust == 1)
datExpr <- datExpr0[keepSamples, ]
nGenes <- ncol(datExpr)
nSamples <- nrow(datExpr)
print(paste("Samples:", nSamples, "Genes(metabolites):", nGenes))

## 2. 选择软阈值功率（构建网络关键步）
powers <- c(1:10, seq(12, 20, 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

# 绘制软阈值结果
par(mfrow = c(1, 2))
cex1 <- 0.9
plot(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (R^2)",
     type = "n",
     main = "Scale independence")
text(sft$fitIndices[, 1],
     -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.9, col = "red")

plot(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = "Mean connectivity")
text(sft$fitIndices[, 1],
     sft$fitIndices[, 5],
     labels = powers, cex = cex1, col = "red")

# 根据图形选择软阈值，比如5
softPower <- 6

## 3. 构建网络并识别模块
# 邻接矩阵和 TOM 相似度
adjacency <- adjacency(datExpr, power = softPower)
TOM <- TOMsimilarity(adjacency)
rownames(TOM) <- colnames(TOM) <- names(datExpr) 
dissTOM <- 1 - TOM

# 层次聚类
geneTree <- hclust(as.dist(dissTOM), method = "average")
plot(geneTree, main = "Metabolite clustering on TOM-based dissimilarity",
     xlab = "", sub = "", labels = FALSE, hang = 0.04)

## 动态剪切模块检测
minModuleSize <- 100  # 模块最小代谢物数，可调整
dynamicMods <- cutreeDynamic(dendro = geneTree,
                             distM = dissTOM,
                             deepSplit = 2,
                             pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors <- labels2colors(dynamicMods)
table(dynamicColors)

# 可视化模块颜色
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Metabolite dendrogram and module colors")

## 4. 模块特征值计算与模块合并
MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
MEs <- MEList$eigengenes
MEDiss <- 1 - cor(MEs)
METree <- hclust(as.dist(MEDiss), method = "average")

plot(METree, main = "Clustering of module eigengenes", xlab = "", sub = "")
MEDissThres <- 0.2  # 模块合并阈值，可调整
abline(h = MEDissThres, col = "red")

merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs

# 可视化合并模块结果
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged modules"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Metabolite dendrogram and merged module colors")

# 最终结果更新
moduleColors <- mergedColors
MEs <- mergedMEs

## 保存模块结果
save(datExpr, moduleColors, MEs, geneTree, TOM,
     file = "CA_WGCNA_metabolite_modules.RData")


## 导出每个模块的网络边与节点
outputDir <- "/Users/liupeihao/Desktop/CTWGCNA"
if (!dir.exists(outputDir)) {dir.create(outputDir, recursive = TRUE)}

threshold <- 0.05
metaboliteNames <- names(datExpr)

uniqueModules <- unique(moduleColors)
for (module in uniqueModules) {
  if (module == "grey") next
  probes <- metaboliteNames[moduleColors == module]
  moduleTOM <- TOM[probes, probes]
  moduleTOM[moduleTOM < threshold] <- 0
  edgeFile <- file.path(outputDir, paste0("CA_CytoscapeInput_edges_", module, ".txt"))
  nodeFile <- file.path(outputDir, paste0("CA_CytoscapeInput_nodes_", module, ".txt"))
  cyt <- exportNetworkToCytoscape(
    moduleTOM,
    edgeFile  = edgeFile,
    nodeFile  = nodeFile,
    weighted  = TRUE,
    threshold = threshold,
    nodeNames = probes,
    altNodeNames = probes,
    nodeAttr  = module  )
  cat("Module", module, "exported to", outputDir, "with", length(probes), "metabolites\n")}

## 代谢物聚类热图数据的绘制与提起
subsetTOM <- TOM
diag(subsetTOM) <- NA

write.csv(subsetTOM, file = "CA_Metabolite_TOM_matrix.csv", quote = FALSE, row.names = TRUE)

module_info <- data.frame(
  Metabolite = colnames(subsetTOM),
  ModuleColor = moduleColors)
write.csv(module_info, "CA_Metabolite_colors.csv", row.names = FALSE)