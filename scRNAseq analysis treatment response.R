# script to analyze and annotate cell clusters in Notch inhibition response of embryos (based on the workflow https://satijalab.org/seurat/archive/v3.1/immune_alignment.html)

library(SingleR)
library(celldex)
library(Seurat)
library(tidyverse)
library(pheatmap)

load("genes.nv2.RObj")

raw.data1 <- Read10X_h5(filename = '2d.DMSO.control/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)
raw.data2 <- Read10X_h5(filename = '2d.LY.treated/filtered_feature_bc_matrix.h5', use.names = TRUE, unique.features = TRUE)

rownames(raw.data1) <- genes$gene_short_name
rownames(raw.data2) <- genes$gene_short_name

# Set up control object
ctrl <- CreateSeuratObject(counts = raw.data1, project = "dmso.48h", min.cells = 5)
ctrl$stim <- "dmso.48h"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)

# Set up stimulated object
stim <- CreateSeuratObject(counts = raw.data2, project = "ly.48h", min.cells = 5)
stim$stim <- "ly.48h"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)

immune.anchors <- FindIntegrationAnchors(object.list = list(ctrl, stim), dims = 1:20)
seurat_obj_dt <- IntegrateData(anchorset = immune.anchors, dims = 1:20)

DefaultAssay(seurat_obj_dt) <- "integrated"

# Run the standard workflow for visualization and clustering
seurat_obj_dt <- ScaleData(seurat_obj_dt, verbose = FALSE)
seurat_obj_dt <- RunPCA(seurat_obj_dt, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
seurat_obj_dt <- RunUMAP(seurat_obj_dt, reduction = "pca", dims = 1:20)
seurat_obj_dt <- FindNeighbors(seurat_obj_dt, reduction = "pca", dims = 1:20)
seurat_obj_dt <- FindClusters(seurat_obj_dt, resolution = 0.6)

# Visualization
p1 <- DimPlot(seurat_obj_dt, reduction = "umap", group.by = "stim")
p2 <- DimPlot(seurat_obj_dt, reduction = "umap", label = TRUE)
plot_grid(p1, p2)

#singleR label transfer
load("Reference.Robj")
reference5 <- as.SingleCellExperiment(reference5)
seurat_obj_dt[["RNA"]] <- JoinLayers(seurat_obj_dt[["RNA"]])
merged_counts <- GetAssayData(seurat_obj_dt, layer = 'counts')
merged_counts <- LayerData(seurat_obj_dt, assay = "RNA",layer = 'counts') #raw counts
merged_counts <- LayerData(seurat_obj_dt, assay = "RNA",layer = 'data') #normalized counts

pred <- SingleR(test = merged_counts,
                ref = reference5,
                labels = reference5$newIDs,de.method = "wilcox")

pred

seurat_obj_dt$singleR.labels <- pred$labels[match(rownames(seurat_obj_dt@meta.data), rownames(pred))]
c7 <- c(`Ectoderm` = "dodgerblue2", `Endoderm` = "green4",`Mesoderm` = "#E31A1C",`Secretory Progenitor` = "#FF7F00",`Cnidocytes` = "black", `Gland Cells` = "gold1",`Neuronal` = "orchid1")
DimPlot(seurat_obj_dt, reduction = 'umap', group.by = 'singleR.labels',cols = c7,pt.size = 1,order = TRUE,split.by = "orig.ident") &NoAxes()
FeaturePlot(object = seurat_obj_dt,features = c("Cadherin1"),pt.size = 1,order = TRUE,reduction = "umap",split.by = "orig.ident") &NoLegend() &NoAxes()