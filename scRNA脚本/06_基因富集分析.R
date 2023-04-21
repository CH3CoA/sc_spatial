#基因富集分析


#单细胞测序如何做GSVA分析1
library(Seurat)
pbmc <- readRDS('pbmc.rds')
DimPlot(pbmc)

#取出B细胞和NK细胞作差异分析
pbmc$celltype <- Idents(pbmc)
pbmc <- subset(pbmc,idents = c('B','NK'))

#准备基因集
#找到B cell相较于NK而言表达量最高的10个基因和表达量最低的10个基因
bcell.marker <- FindMarkers(pbmc,ident.1 = 'B',ident.2 = 'NK')
library(dplyr)
top10gene <- bcell.marker %>% top_n(n = 10, wt = avg_log2FC) %>% rownames()
bottom10gene <- bcell.marker %>% top_n(n = -10, wt = avg_log2FC) %>% rownames()

#查看高表达的10个基因和低表达的10个基因
top10gene
bottom10gene

#构建基因集用于GSVA
mygeneset <- cbind(top10gene,bottom10gene)
mygeneset <- as.data.frame(mygeneset)

#表达量，你得这么取,取出表达矩阵
gene.expr <-  as.matrix(pbmc[["RNA"]]@data)
dim(gene.expr)

#GSVA 分析反而是最容易的
#如果你这里看不明白，参考一下我们前面的几篇推送
library(GSVA)
gsva.result <- GSVA::gsva(gene.expr, mygeneset,kcdf='Gaussian')

#处理一下表达矩阵
mymatrix <- t(gsva.result)
mymatrix <- cbind(mymatrix,as.data.frame(pbmc$celltype))
colnames(mymatrix)[ncol(mymatrix)] <- 'celltype'
head(mymatrix)

#把这两个通路的GSVA结果画成箱线图
library(ggplot2)
## Warning: package 'ggplot2' was built under R version 4.0.5
mytop <- ggplot2::ggplot(mymatrix,aes(x=celltype,y=top10gene,fill=celltype))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  labs(title ='B Cell Top10gene')


mybottom <- ggplot2::ggplot(mymatrix,aes(x=celltype,y=bottom10gene,fill=celltype))+
  geom_boxplot(alpha=0.7)+
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name="Celltype")+
  labs(title ='B cell bottom10gene')
library(patchwork)
mytop|mybottom









#富集分析 clusterprofiler
suppressWarnings(suppressMessages(if(!require(rvcheck))devtools::install_version("rvcheck", version = "0.1.8", repos = "http://cran.us.r-project.org")))
suppressWarnings(suppressMessages(if(!require(clusterProfiler))BiocManager::install("clusterProfiler")))
suppressWarnings(suppressMessages(if(!require(org.Mm.eg.db))BiocManager::install("org.Mm.eg.db")))
suppressWarnings(suppressMessages(if(!require(org.Hs.eg.db))BiocManager::install("org.Hs.eg.db")))
suppressWarnings(suppressMessages(library(dplyr)))

#准备基因和绘图函数
bcell.marker <- read.csv('bcell.marker.CSV',row.names = 1)
head(bcell.marker)

#取出100个p_val_adj最小的基因
mygene <- bcell.marker %>% top_n(n = -100, wt = p_val_adj) %>% rownames()
mygene[1:10]

#这里记录一下clusterprofile的ID转换功能
gene.df <- bitr(mygene,fromType="SYMBOL",toType=c("ENTREZID","ENSEMBL"),
                OrgDb = org.Hs.eg.db)# org.Mm.eg.db小鼠

#查看gene.df
gene.df

#可视化的函数（最后运行，绘图）
erich2plot <- function(data4plot){
  library(ggplot2)
  data4plot <- data4plot[order(data4plot$qvalue,decreasing = F)[1:30],]
  data4plot$BgRatio<-
    apply(data4plot,1,function(x){
      as.numeric(strsplit(x[3],'/')[[1]][1])
    })/apply(data4plot,1,function(x){
      as.numeric(strsplit(x[4],'/')[[1]][1])
    })
  
  p <- ggplot(data4plot,aes(BgRatio,Description))
  p<-p + geom_point()
  
  pbubble <- p + geom_point(aes(size=Count,color=-1*log10(qvalue)))
  
  pr <- pbubble + scale_colour_gradient(low="#90EE90",high="red") + 
    labs(color=expression(-log[10](qvalue)),size="observed.gene.count", 
         x="Richfactor", y="term.description",title="Enrichment Process")
  
  pr <- pr + theme_bw()
  pr
}

erich2plot(ekegg@result)

#保存为pdf
pdf(file = 'kegg.pdf',width = 10,height = 10)
erich2plot(ekegg@result)
dev.off()

#KEGG 与 GO 富集分析
#使用在线数据,可能会受域名限制（人has，小鼠mmu）
ekegg <- enrichKEGG(unique(gene.df$ENTREZID), organism='hsa',
                    pvalueCutoff=0.05,pAdjustMethod='BH',qvalueCutoff=0.2,
                    minGSSize=10,maxGSSize=500,use_internal_data=F)

#分析结果在list的result中（取出查看：富集分析的ID、Description、GeneRatio、BgRatio、p值等）
#查看ekegg@result
lalala<-ekegg@result
#将ENTREZID ID转换成SYMBOL
ekegg <- setReadable(ekegg,'org.Hs.eg.db','ENTREZID')
write.csv(ekegg@result,'kegg.result.csv')

#查看ekegg@result
lalala<-ekegg@result



#######GO over-representation test###########
ggoMF <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "MF",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoMF@result,'GO.MF.result.csv')

pdf(file = 'GO.MF.pdf',width = 10,height = 10)
erich2plot(ggoMF@result)
dev.off()


ggoCC <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "CC",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoCC@result,'GO.CC.result.csv')

pdf(file = 'ggoCC.pdf',width = 10,height = 10)
erich2plot(ggoCC@result)
dev.off()


ggoBP <- enrichGO(gene.df$SYMBOL, org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                  pvalueCutoff = 0.05, pAdjustMethod = "BH", qvalueCutoff = 0.2, minGSSize = 10,
                  maxGSSize = 500, readable = FALSE, pool = FALSE)
write.csv(ggoBP@result,'GO.BP.result.csv')

pdf(file = 'ggoBP.pdf',width = 10,height = 10)
erich2plot(ggoBP@result)
dev.off()



#GSEA 分析
bcell.marker$SYMBOL <- rownames(bcell.marker)
genetable <- merge(gene.df,bcell.marker,by='SYMBOL')  
geneList <- genetable$avg_log2FC
names(geneList) <- genetable$ENTREZID

geneList <- sort(geneList,decreasing = T)
egseGO <- gseGO(geneList,ont = "MF", org.Hs.eg.db, exponent = 1, nPerm = 1000, 
                minGSSize = 10, maxGSSize = 500, pvalueCutoff = 1, pAdjustMethod = "BH", 
                verbose = TRUE, seed = FALSE, by = "fgsea")#换成CC   BP

kk2 <- gseKEGG(geneList=geneList,
               organism = "hsa",#hsa => human   , mmu => mouse #可用物种列表：   
               #https://www.genome.jp/kegg/catalog/org_list.html
               nPerm=1000,
               minGSSize = 120,
               pvalueCutoff = 1,
               verbose = FALSE)

for (i in 1:length(egseGO$Description)) {
  plot1 <- gseaplot(egseGO,geneSetID = i,title = egseGO$Description[i])
  ggsave( paste0(egseGO$Description[i],".pdf"),plot = plot1,height = 10,width = 10)
}

plot1


#clusterProfiler 包里的一些默认作图方法，例如
barplot(ekegg,)  #富集柱形图
dotplot(ekegg, showCategory=20)  #富集气泡图
cnetplot(ekegg) #网络图展示富集功能和基因的包含关系
emapplot(ekegg) #网络图展示各富集功能之间共有基因关系
heatplot(ekegg) #热图展示富集功能和基因的包含关系
