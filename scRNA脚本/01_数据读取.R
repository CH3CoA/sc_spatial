#数据读取
#数据下载
#https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi

library(multtest)
if(!require(multtest))install.packages("multtest")
if(!require(Seurat))install.packages("Seurat")
if(!require(dplyr))install.packages("dplyr")
if(!require(mindr))install.packages("mindr")
if(!require(tidyverse))install.packages("tidyverse")

#####自动读取cellranger(LINUX)输出的feature barcode matric

#清空
rm(list = ls())
setwd('E:/Python/Spatial_Transcriptomics')
pbmc.data <- Read10X(data.dir = "data/GSE181276_genes.counts_for_GEO_uploading.txt/")
#自动读取10X的数据，是一些tsv与mtx文件
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "biomamba")

#设置参数，可忽略（基因至少在3个细胞内表达，细胞至少表达200个基因，才会被保留）
pbmc<CreateSeuratObject(counts=pbmc.data,project="biomamba",min.cells=3,min.features=200)


####仅有一个稀疏矩阵时的读取方法#####
matrix_data <- read.table("data/GSE181276_genes.counts_for_GEO_uploading.txt/GSE181276_genes.counts_for_GEO_uploading.txt", sep="\t", header=T, row.names=1)

dim(matrix_data)
#13714个基因 2700个细胞

#Seurat包
seurat_obj <- CreateSeuratObject(counts = matrix_data)

######读取RDS文件############
rm(list = ls())
pbmc <- readRDS("panc8.rds")
saveRDS(pbmc,"pbmc.rds")


#############
library(tidyverse)
str(pbmc)
library(mindr)
(out <- capture.output(str(pbmc)))
out2 <- paste(out, collapse="\n")
mm(gsub("\\.\\.@","# ",gsub("\\.\\. ","#",out2)),type ="text",root= "Seurat")
