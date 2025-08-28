library(Seurat)
library(SingleR)
library(SingleCellExperiment)

#singleR label transfer for gastruloids atlas 
load('object from script scRNAseq analysis gastruloids')
load("Reference.Robj")
reference5 <- as.SingleCellExperiment(reference5)
merged_counts <- GetAssayData(seurat_obj_g, layer = 'counts')
merged_counts <- LayerData(seurat_obj_g, assay = "RNA",layer = 'counts') #raw counts
merged_counts <- LayerData(seurat_obj_g, assay = "RNA",layer = 'data') #normalized counts

pred <- SingleR(test = merged_counts,
                ref = reference5,
                labels = reference5$newIDs,de.method = "wilcox")

pred

seurat_obj_g$singleR.labels <- pred$labels[match(rownames(seurat_obj_g@meta.data), rownames(pred))]
DimPlot(seurat_obj_g, reduction = 'umap', group.by = 'singleR.labels',cols = c7,pt.size = 1,order = TRUE,split.by = "newIdents") &NoAxes()

#singleR label transfer for drug treatment response
load('object from script scRNAseq analysis treatment response')
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