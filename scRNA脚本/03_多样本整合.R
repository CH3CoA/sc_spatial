#多样本整合
###########单纯的merge#################
library(Seurat)
library(multtest)
library(dplyr)
library(ggplot2)
library(patchwork)

##########准备用于拆分的数据集##########
#pbmc <- subset(pbmc, downsample = 50)每个cluster随机取50个细胞，小的数据集
ifnb <- readRDS('pbmcrenamed.rds')
ifnb.list <- SplitObject(ifnb, split.by = "group")

unique(ifnb$group)#查看ifnb里面iroup的变量

C57 <- ifnb.list$C57
AS1 <- ifnb.list$AS1

######简单merge######## 
#不具有去批次效应功能
pbmc <- merge(C57, y = c(AS1), add.cell.ids = c("C57", "AS1"), project = "ALL")
pbmc

head(colnames(pbmc))
#组别信息

unique(sapply(X = strsplit(colnames(pbmc), split = "_"), FUN = "[", 1))
#还原组别

table(pbmc$orig.ident)
统计细胞数量

##############anchor###############另一种方法
library(Seurat)
library(tidyverse)
### testA ----
myfunction1 <- function(testA.seu){
  testA.seu <- NormalizeData(testA.seu, normalization.method = "LogNormalize", scale.factor = 10000)
  testA.seu <- FindVariableFeatures(testA.seu, selection.method = "vst", nfeatures = 2000)
  return(testA.seu)
}
C57 <- myfunction1(C57)
AS1 <- myfunction1(AS1)

### Integration ----
testAB.anchors <- FindIntegrationAnchors(object.list = list(C57,AS1), dims = 1:20)
testAB.integrated <- IntegrateData(anchorset = testAB.anchors, dims = 1:20)

#需要注意的是：上面的整合步骤相对于harmony整合方法，对于较大的数据集（几万个细胞）
#非常消耗内存和时间，大约9G的数据32G的内存就已经无法运行；
#当存在某一个Seurat对象细胞数很少（印象中200以下这样子），
#会报错，这时建议用第二种整合方法

DefaultAssay(testAB.integrated) <- "integrated"

# # Run the standard workflow for visualization and clustering
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.5)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:30)
p1<- DimPlot(testAB.integrated,label = T,split.by = 'group')#integrated

DefaultAssay(testAB.integrated) <- "RNA"
testAB.integrated <- ScaleData(testAB.integrated, features = rownames(testAB.integrated))
testAB.integrated <- RunPCA(testAB.integrated, npcs = 50, verbose = FALSE)
testAB.integrated <- FindNeighbors(testAB.integrated, dims = 1:30)
testAB.integrated <- FindClusters(testAB.integrated, resolution = 0.5)
testAB.integrated <- RunUMAP(testAB.integrated, dims = 1:30)
testAB.integrated <- RunTSNE(testAB.integrated, dims = 1:30)

p2 <- DimPlot(testAB.integrated,label = T,split.by = 'group')

p1|p2#对比p1和p2

###########harmony 速度快、内存少################
library(devtools)
if(!require(harmony))devtools::install_github("immunogenomics/harmony")
test.seu <- pbmc
test.seu <-  test.seu%>%
  Seurat::NormalizeData() %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>% 
  ScaleData()
test.seu <- RunPCA(test.seu, npcs = 50, verbose = FALSE)

#####run 到PCA再进行harmony，相当于降维########
test.seu=test.seu %>% RunHarmony("group", plot_convergence = TRUE)

test.seu <- test.seu %>% 
  RunUMAP(reduction = "harmony", dims = 1:30) %>% 
  FindNeighbors(reduction = "harmony", dims = 1:30) %>% 
  FindClusters(resolution = 0.5) %>% 
  identity()

test.seu <- test.seu %>% 
  RunTSNE(reduction = "harmony", dims = 1:30)

p3 <- DimPlot(test.seu, reduction = "tsne", group.by = "group", pt.size=0.5)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)
p4 <- DimPlot(test.seu, reduction = "tsne", group.by = "ident",   pt.size=0.5, label = TRUE,repel = TRUE)+theme(
  axis.line = element_blank(),
  axis.ticks = element_blank(),axis.text = element_blank()
)

p3|p4
