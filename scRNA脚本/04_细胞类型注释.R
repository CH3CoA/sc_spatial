#细胞类型注释
#准备工作，先进行预处理，作到细胞注释前的步骤
if(T){rm(list = ls())
  
  if (!require("Seurat"))install.packages("Seurat")
  if (!require("BiocManager", quietly = TRUE))install.packages("BiocManager")
  if (!require("multtest", quietly = TRUE))install.packages("multtest")
  if (!require("dplyr", quietly = TRUE))install.packages("dplyr")
  download.file('https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz','pbmc3k_filtered_gene_bc_matrices.tar.gz')
  library(R.utils)
  gunzip('pbmc3k_filtered_gene_bc_matrices.tar.gz')
  untar('pbmc3k_filtered_gene_bc_matrices.tar')
  pbmc.data <- Read10X('filtered_gene_bc_matrices/hg19/') 
  pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
  pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)      
  pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
  top10 <- head(VariableFeatures(pbmc), 10) 
  pbmc <- ScaleData(pbmc, features =  rownames(pbmc)) 
  pbmc <- ScaleData(pbmc, vars.to.regress = "percent.mt") 
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  print(pbmc[["pca"]], dims = 1:5, nfeatures = 5) 
  VizDimLoadings(pbmc, dims = 1:2, reduction = "pca") 
  DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE) 
  DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE) 
  pbmc <- JackStraw(pbmc, num.replicate = 100) 
  pbmc <- ScoreJackStraw(pbmc, dims = 1:20) 
  JackStrawPlot(pbmc, dims = 1:15)
  pbmc <- FindNeighbors(pbmc, dims = 1:10)
  pbmc <- FindClusters(pbmc, resolution = 0.5) 
  head(Idents(pbmc), 5) 
  pbmc <- RunUMAP(pbmc, dims = 1:10) 
  pbmc <- RunTSNE(pbmc, dims = 1:10) 
  DimPlot(pbmc, reduction = "umap", label = TRUE)
  pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) 
  library(dplyr)
  pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC) 
}


#方法一：查数据
library(dplyr)
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
#热图
DoHeatmap(pbmc, features = top10$gene) + NoLegend()#展示前10个标记基因的热图
#小提琴图
VlnPlot(pbmc, features = top10$gene[1:20],pt.size=0)
#VlnPlot(pbmc, features = top10$gene[1:20])
#单细胞聚类图
DimPlot(pbmc,label = T)

#通过标记基因及文献，可以人工确定各分类群的细胞类型，则可以如下手动添加细胞群名称
bfreaname.pbmc <- pbmc
new.cluster.ids <- c("Naive CD4 T", "Memory CD4 T", "CD14+ Mono", "B", "CD8 T", "FCGR3A+ Mono", 
                     "NK", "DC", "Platelet")
#帮助单细胞测序进行注释的数据库：
#http://mp.weixin.qq.com/s?__biz=MzI5MTcwNjA4NQ==&mid=2247502903&idx=2&sn=fd21e6e111f57a4a2b6c987e391068fd&chksm=ec0e09bddb7980abf038f62d03d3beea6249753c8fba69b69f399de9854fc208ca863ca5bc23&mpshare=1&scene=24&srcid=1110SJhxDL8hmNB5BThrgOS9&sharer_sharetime=1604979334616&sharer_shareid=853c5fb0f1636baa0a65973e8b5db684#rd
#cellmarker: http://biocc.hrbmu.edu.cn/CellMarker/index.jsp
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

#方法二：通过singleR进行注释 #### #这一方法很鸡肋，主要是由于很难找到合适的参考数据，但是对于免疫细胞的注释还是有可观效果的
if(!require(SingleR))BiocManager::install("SingleR")

if(!require(matrixStats))BiocManager::install('matrixStats')

if(!require(celldex))BiocManager::install('celldex')#有一些版本冲突，需要重新安装一些包

cg<-ImmGenData()#选取我们要使用的免疫细胞参考数据集

assay_for_SingleR <- GetAssayData(bfreaname.pbmc, slot="data")#取出样本中的表达序列
predictions <- SingleR(test=assay_for_SingleR, 
                       ref=cg, labels=cg$label.main)
#以kidney中提取的阵列为输入数据，以小鼠的阵列作为参考，predict细胞类型

table(predictions$labels)#看看都注释到了哪些细胞

cellType=data.frame(seurat=bfreaname.pbmc@meta.data$seurat_clusters,
                    predict=predictions$labels)#得到seurat中编号与预测标签之间的关系
sort(table(cellType[,1]))

table(cellType[,1:2])#访问celltyple的2~3列

#可以看出 singleR如果没有合适的数据集，得到的结果有多拉跨了吧
#这里就没必要再往下重命名了

当你有合适的参考数据集时可以用方法三方法四进行注释
#方法三：自定义singleR的注释
#我们走个极端，拿它自己作为自己的参考数据集，看看注释的准不准
##########利用singleR构建自己的数据作为参考数据集########
library(SingleR)
library(Seurat)
library(ggplot2)
if(!require(textshape))install.packages('textshape')
if(!require(scater))BiocManager::install('scater')
if(!require(SingleCellExperiment))BiocManager::install('SingleCellExperiment')
library(dplyr)

# 读入scRNA数据 -------
myref <- pbmc##这里为了检测，我们将参考数据集与目标数据集用同一个数据进行测试
myref$celltype <- Idents(myref)
table(Idents(myref))

# 读入参考数据集 -------
Refassay <- log1p(AverageExpression(myref, verbose = FALSE)$RNA)#求
#Ref <- textshape::column_to_rownames(Ref, loc = 1)#另一种得到参考矩阵的办法
head(Refassay)#看看表达矩阵长啥样

ref_sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=Refassay))
#参考数据集需要构建成一个SingleCellExperiment对象
ref_sce=scater::logNormCounts(ref_sce)

logcounts(ref_sce)[1:4,1:4]

colData(ref_sce)$Type=colnames(Refassay)
ref_sce#构建完成

testdata <- GetAssayData(bfreaname.pbmc, slot="data")
pred <- SingleR(test=testdata, ref=ref_sce, 
                labels=ref_sce$Type,
                #clusters = scRNA@active.ident
)
table(pred$labels)

head(pred)

as.data.frame(table(pred$labels))

#pred@listData[["scores"]] #预测评分，想看看结构的可以自己看看
#同上，我们找一下seurat中类群与注释结果直接的关系
cellType=data.frame(seurat=bfreaname.pbmc@meta.data$seurat_clusters,
                    predict=pred$labels)#得到seurat中编号与预测标签之间的关系
sort(table(cellType[,1]))

table(cellType[,1:2])#访问celltyple的2~3列

lalala <- as.data.frame(table(cellType[,1:2]))
finalmap <- lalala %>% group_by(seurat) %>% top_n(n = 1, wt = Freq)#找出每种seurat_cluster注释比例最高的对应类型
finalmap <-finalmap[order(finalmap$seurat),]$predict#找到seurat中0：8的对应预测细胞类型
print(finalmap)

testname <- bfreaname.pbmc
new.cluster.ids <- as.character(finalmap)
names(new.cluster.ids) <- levels(testname)
testname <- RenameIdents(testname, new.cluster.ids)

p1 <- DimPlot(pbmc,label = T)
p2 <- DimPlot(testname,label = T)#比较一下测试数据与参考数据集之间有没有偏差
p1|p2#完美，无差别注释，当然了，我们这个参考数据用的比较极端





######利用seurat内置的原先用于细胞整合的功能，将参考数据与待注释数据进行映射处理
library(Seurat)
pancreas.query <- bfreaname.pbmc#待注释数据
pancreas.anchors <- FindTransferAnchors(reference = pbmc, query = pancreas.query,
                                        dims = 1:30)

pbmc$celltype <- Idents(pbmc)
predictions <- TransferData(anchorset = pancreas.anchors, refdata = pbmc$celltype,
                            dims = 1:30)

pancreas.query <- AddMetaData(pancreas.query, metadata = predictions)
#把注释加回原来的数据集
pancreas.query$prediction.match <- pancreas.query$predicted.id
table(pancreas.query$prediction.match)

Idents(pancreas.query)<- 'prediction.match'

p1 <- DimPlot(pbmc,label = T)
p2 <- DimPlot(pancreas.query,label = T)#比较一下测试数据与参考数据集之间有没有偏差
p1|p2#完美，无差别注释，当然了，我们这个参考数据用的比较极端





#### 最近发现代码运行的最大障碍是各个packages的版本冲突问题，这里列出本次分析环境中的所有信息，报错时可以参考
sessionInfo()
