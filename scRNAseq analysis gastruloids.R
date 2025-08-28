#script for scRNAseq analysis of gastruloid atlas 
library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

load("genes.nv2.RObj")

raw.data1 <- Read10X(data.dir = 'aggregates/24h.gast(ACME)/')
raw.data1=raw.data1[c(-24540,-24544),]
raw.data2 <- Read10X(data.dir = 'aggregates/4h/')
raw.data3 <- Read10X(data.dir = 'aggregates/12h/')
raw.data4 <- Read10X(data.dir = 'aggregates/24h/')
raw.data5 <- Read10X_h5(filename = 'aggregates/24h_rep/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
raw.data6 <- Read10X(data.dir = 'aggregates/48h/')
raw.data7 <- Read10X_h5(filename = 'aggregates/48h_rep/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

rownames(raw.data1) <- genes$gene_short_name
rownames(raw.data2) <- genes$gene_short_name
rownames(raw.data3) <- genes$gene_short_name
rownames(raw.data4) <- genes$gene_short_name
rownames(raw.data5) <- genes$gene_short_name
rownames(raw.data6) <- genes$gene_short_name
rownames(raw.data7) <- genes$gene_short_name

wt.24h  <- CreateSeuratObject(counts = raw.data1, project = "wt.24h") 
agg.4h  <- CreateSeuratObject(counts = raw.data2, project = "agg.4h") 
agg.12h  <- CreateSeuratObject(counts = raw.data3, project = "agg.12h") 
agg.24h_1 <- CreateSeuratObject(counts = raw.data4, project = "agg.24h_1") 
agg.24h_2 <- CreateSeuratObject(counts = raw.data5, project = "agg.24h_2") 
agg.48h_1  <- CreateSeuratObject(counts = raw.data6, project = "agg.48h_1") 
agg.48h_2  <- CreateSeuratObject(counts = raw.data7, project = "agg.48h_2") 

agg.24h <- merge(agg.24h_1,agg.24h_2)
agg.48h <- merge(agg.48h_1,agg.48h_2)

levels(wt.24h@meta.data$orig.ident) <- 'wt.24h'
levels(agg.4h@meta.data$orig.ident) <- 'agg.4h'
levels(agg.12h@meta.data$orig.ident) <- 'agg.12h'
levels(agg.24h@meta.data$orig.ident) <- 'agg.24h'
levels(agg.48h@meta.data$orig.ident) <- 'agg.48h'

VlnPlot(wt.24h, features = c('nFeature_RNA','nCount_RNA'))
VlnPlot(agg.4h, features = c('nFeature_RNA','nCount_RNA'))
VlnPlot(agg.12h, features = c('nFeature_RNA','nCount_RNA'))
VlnPlot(agg.24h, features = c('nFeature_RNA','nCount_RNA'))
VlnPlot(agg.48h, features = c('nFeature_RNA','nCount_RNA'))

wt.24h <- subset(x = wt.24h, subset = nFeature_RNA > 500 & nCount_RNA < 15000)
agg.4h <- subset(x = agg.4h, subset = nFeature_RNA > 500 & nCount_RNA < 15000)
agg.12h <- subset(x = agg.12h, subset = nFeature_RNA > 500 & nCount_RNA < 15000)
agg.24h <- subset(x = agg.24h, subset = nFeature_RNA > 500 & nCount_RNA < 15000)
agg.48h <- subset(x = agg.48h, subset = nFeature_RNA > 500 & nCount_RNA < 15000)

wt.24h <- RenameCells(wt.24h, add.cell.id = "wt.24h")
agg.4h <- RenameCells(agg.4h, add.cell.id = "agg.4h")
agg.12h <- RenameCells(agg.12h, add.cell.id = "agg.12h")
agg.24h <- RenameCells(agg.24h, add.cell.id = "agg.24h")
agg.48h <- RenameCells(agg.48h, add.cell.id = "agg.48h")

seurat_obj_g <- merge(x = wt.24h, y = list(agg.4h,agg.12h,agg.24h,agg.48h))
seurat_obj_g <- NormalizeData(seurat_obj_g)
seurat_obj_g <- FindVariableFeatures(seurat_obj_g)
seurat_obj_g <- ScaleData(seurat_obj_g)
seurat_obj_g <- RunPCA(seurat_obj_g)
seurat_obj_g <- IntegrateLayers(object = seurat_obj_g, method = HarmonyIntegration, orig.reduction = "pca",new.reduction = "harmony")
# now that integration is complete, rejoin layers
seurat_obj_g[["RNA"]] <- JoinLayers(seurat_obj_g[["RNA"]])
seurat_obj_g <- FindNeighbors(seurat_obj_g, reduction = "harmony", dims = 1:30)
seurat_obj_g <- FindClusters(seurat_obj_g, resolution = 0.4)
seurat_obj_g <- RunUMAP(seurat_obj_g, dims = 1:30, reduction = "harmony")
seurat_obj_g$orig.ident <- factor(seurat_obj_g$orig.ident,levels=c("wt.24h","agg.4h","agg.12h","agg.24h_1","agg.24h_2"
                                                               ,"agg.48h_1","agg.48h_2"))
Idents(seurat_obj_g) <- seurat_obj_g$orig.ident
seurat_obj_g <- RenameIdents(object = seurat_obj_g,`wt.24h` = "wt.24h",`agg.4h` = "agg.4h",`agg.12h` = "agg.12h",`agg.24h_1` = "agg.24h",`agg.24h_2` = "agg.24h",`agg.48h_1` = "agg.48h",`agg.48h_2` = "agg.48h")
seurat_obj_g$newIdents <- Idents(seurat_obj_g)
Idents(seurat_obj_g) <- 'seurat_clusters'
c20 <- c(`0` = "dodgerblue2", `1` = "#E31A1C",`2` = "green4",`3` = "#6A3D9A",`4` = "#FF7F00",`5` = "black", `6` = "gold1",`7` = "skyblue2",`8` = "palegreen2",`9` = "gray70", `10` = "khaki2",`11` = "maroon", `12` = "orchid1", `13` = "deeppink1",`14` =  "blue1", `15` = "steelblue4",`16` = "darkturquoise",`17` =  "green1", `18` = "yellow4", `19` = "brown")
DimPlot(seurat_obj_g,pt.size = 1,split.by = "newIdents",reduction = "umap",order = TRUE,cols = c20) &NoAxes() 
FeaturePlot(seurat_obj_g,features = c("GAPR1-like-7"),pt.size = 1,order = TRUE,split.by = "newIdents",reduction = "umap") &NoAxes() &NoLegend()