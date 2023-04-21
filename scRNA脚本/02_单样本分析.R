#单样本分析
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(patchwork))install.packages("patchwork")
if(!require(R.utils))install.packages("R.utils")

#下载并解压数据
download.file('https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz','my.pbmc.gz')
untar(gunzip("my.pbmc.gz"))

#读入数据并创建Seurat分析对象
pbmc.data <- Read10X(data.dir = "filtered_gene_bc_matrices/hg19/")
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc
ncol(pbmc)
ncol(pbmc.data)
lalala <- as.data.frame(pbmc[["RNA"]]@counts)
write.table(lalala,'mycount.txt',sep = '\t')#表达矩阵可以这么存出来

#质控数据及可视化（feature画的纵坐标，ncol每一行画几张图；mt¬_线粒体RNA的比例,过高_细胞活性不好或细胞线粒体含量较高，过低_）
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-") 
#鼠源的需要换成mt
#人源MT
head(pbmc@meta.data,5)
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
if(!require(patchwork))install.packages("patchwork")
CombinePlots(plots = list(plot1, plot2))
VlnPlot(pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)


#过滤参数（质控：过少_细胞活性不好或检测深度不够，过多_双细胞或多细胞）
pbmc <- subset(pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)   
ncol(as.data.frame(pbmc[["RNA"]]@counts))
#查看剩余细胞数
pbmc

#计算、分群（Normalize每一个细胞标准化）
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
#寻找高变基因
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(pbmc), 10)
# Identify the 10 most highly variable genes
#高变基因可视化
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1 + plot2

#scaleData每一个基因在细胞中的标准化，防止细胞分群中某一基因表达离群造成的差异
pbmc <- ScaleData(pbmc, features = rownames(pbmc))

## Centering and scaling data matrix
#pbmc <- ScaleData(pbmc) ##only4VariableFeatures
#默认高变基因进行scale，如果是全部数据scale选取下面一个进行scale

#PCA分析，初步降维
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
#PCA结果
print(pbmc[["pca"]], dims = 1:5, nfeatures = 5)
#点图形式展示PCA
VizDimLoadings(pbmc, dims = 1:2, reduction = "pca")
#投影的降维图
DimPlot(pbmc, reduction = "pca")

#可删除
DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)
#前15个维度的PCA表现
DimHeatmap(pbmc, dims = 1:15, cells = 500, balanced = TRUE)

#高阶PCA分析及可视化
pbmc <- JackStraw(pbmc, num.replicate = 100)#时间较长
pbmc <- ScoreJackStraw(pbmc, dims = 1:20)
JackStrawPlot(pbmc, dims = 1:20)#可视化

#可视化
ElbowPlot(pbmc)

#选取前十个PCA
pbmc <- FindNeighbors(pbmc, dims = 1:10)
#细胞分群
pbmc <- FindClusters(pbmc, resolution = 0.5)
# resolution在0.4至1.2之间，resolution越大所获得的细胞类群越多，根据测序数据实时修改，反复调试

#分群后的可视化：tsne+umap
pbmc <- RunUMAP(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "umap")
#umap适合大型数据

pbmc <- RunTSNE(pbmc, dims = 1:10)
DimPlot(pbmc, reduction = "tsne")
#tsne小型数据，但分离效果比umap好

#寻找marker基因并对cluster进行重命名
cluster5.markers <- FindMarkers(pbmc, ident.1 = 5, ident.2 = c(0, 3), min.pct = 0.25)#一般不做
#cluster5细胞类型的marker基因

pbmc.markers <- FindAllMarkers(pbmc, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#每个细胞类型的marker基因

#高表达marker基因作为初始数值传递给管道符，选出表达量前几的基因
if(!require(dplyr))install.packages("dplyr")
pbmc.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)

VlnPlot(pbmc, features = c("NKG7", "PF4"), slot = "counts", log = TRUE)

FeaturePlot(pbmc, features = c("MS4A1", "GNLY", "CD3E", "CD14", "FCER1A", "FCGR3A", "LYZ"))

#热图展示（找每一群细胞高表达marker前十）
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(pbmc, features = top10$gene) + NoLegend()

#http://biocc.hrbmu.edu.cn/CellMarker/（通过数据可或文献根据找marker对应的细胞）

#导出top10的csv格式
write.csv(top10,'top10.csv')
# GSE81076_D2_3_7_10_17.txt
new.cluster.ids <- c("Acinar cell", " Meiotic prophase fetal germ cell or SLC16A7+cell", " Epithelial cell or Pancreatic ductal stem cell", " Alpha cell", " Beta cell", " Acinar cell", "NA", " Delta cell", " Granulocyte-monocyte progenitor or Mast cell ")






new.cluster.ids <- c("Naive CD4 T", "CD14+ Mono", "Memory CD4 T", "B", "CD8 T", "FCGR3A+ Mono", "NK", "DC", "Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc, new.cluster.ids)
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()

sessionInfo()

#保存为RDS文件
saveRDS(pbmc,'pbmc.rds')
pbmc<- readRDS('pbmc.rds')
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
