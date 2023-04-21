#组间差异分析及可视化
#读入并检查数据
library(Seurat)

library(dplyr)

pbmc <- readRDS('pbmcrenamed.rds')
pbmc

DimPlot(pbmc)

names(pbmc@meta.data)

unique(pbmc$group)

DimPlot(pbmc,split.by = 'group')

#差异分析
?Findermarker
Idents(pbmc)
#不用运行，查看目前Idents，目前Idents是细胞类型
pbmc$celltype.group <- paste(pbmc$celltype, pbmc$group, sep = "_")
pbmc$celltype <- Idents(pbmc)
Idents(pbmc) <- "celltype.group"

mydeg <- FindMarkers(pbmc,ident.1 = 'VSMC_AS1',ident.2 = 'VSMC_C57', verbose = FALSE, test.use = 'wilcox',min.pct = 0.1)
head(mydeg)

#解放生产力 通过循环自动计算差异基因
cellfordeg<-levels(pbmc$celltype)
for(i in 1:length(cellfordeg)){
  CELLDEG <- FindMarkers(pbmc, ident.1 = paste0(cellfordeg[i],"_P3"), ident.2 = paste0(cellfordeg[i],"_AS1"), verbose = FALSE)
  write.csv(CELLDEG,paste0(cellfordeg[i],".CSV"))
}
list.files()

#可视化方法
library(dplyr)
top10 <- CELLDEG  %>% top_n(n = 10, wt = avg_log2FC) %>% row.names()
top10

pbmc <- ScaleData(pbmc, features =  rownames(pbmc))

DoHeatmap(pbmc,features = top10,size=3)

Idents(pbmc) <- "celltype"
VlnPlot(pbmc,features = top10,split.by = 'group',idents = 'EC')

FeaturePlot(pbmc,features = top10,split.by = 'group')

#DotPlot(pbmc,features = top10,split.by ='group')#默认只有两种颜色
DotPlot(pbmc,features = top10,split.by ='group',cols = c('blue','yellow','pink'))

#提取表达量，用ggplot2 DIY一个箱线图
####提取表达量#######
mymatrix <- as.data.frame(pbmc@assays$RNA@data)
mymatrix2<-t(mymatrix)%>%as.data.frame()
mymatrix2[,1]<-pbmc$celltype
colnames(mymatrix2)[1] <- "celltype"

mymatrix2[,ncol(mymatrix2)+1]<-pbmc$group
colnames(mymatrix2)[ncol(mymatrix2)] <- "group"

#绘图
library(ggplot2)
p1<- ggplot2::ggplot(mymatrix2,aes(x=celltype,y=Thbs1,fill=group))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange','pink'))
p1

#另一种提取方法
Idents(pbmc) <- colnames(pbmc)
mymatrix <- log1p(AverageExpression(pbmc, verbose = FALSE)$RNA)
mymatrix2<-t(mymatrix)%>%as.data.frame()
mymatrix2[,1]<-pbmc$celltype
colnames(mymatrix2)[1] <- "celltype"

mymatrix2[,ncol(mymatrix2)+1]<-pbmc$group
colnames(mymatrix2)[ncol(mymatrix2)] <- "group"

library(ggplot2)
p2<- ggplot2::ggplot(mymatrix2,aes(x=celltype,y=Thbs1,fill=group))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  scale_fill_manual(values = c('DeepSkyBlue','Orange','pink'))
p2

p1|p2

