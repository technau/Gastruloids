library(hdWGCNA)
library(WGCNA)
library(Seurat)
library(patchwork)

load("genes.nv2.RObj")
load('object from scRNAseq analysis treatment response AFTER singleR.label transfer and annotation')
Idents(seurat_obj_dt) <- "seurat_clusters"

seurat_obj_dt <- SetupForWGCNA(
  seurat_obj_dt,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included
  wgcna_name = "tutorial" # the name of the hdWGCNA experiment
)

seurat_obj_dt <- MetacellsByGroups(
  seurat_obj = seurat_obj_dt,
  group.by = c("seurat_clusters", "orig.ident"), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'umap', # select the dimensionality reduction to perform KNN on
  k = 25, # nearest-neighbors parameter
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'seurat_clusters' # set the Idents of the metacell seurat object
)

seurat_obj_dt <- NormalizeMetacells(seurat_obj_dt)
seurat_obj_dt <- ScaleMetacells(seurat_obj = seurat_obj_dt)

seurat_obj_dt <- SetDatExpr(
  seurat_obj_dt,
  group_name = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13"), # the name of the group of interest in the group.by column
  group.by='seurat_clusters', # the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', # using RNA assay
  layer = 'data' # using normalized data
)

seurat_obj_dt <- TestSoftPowers(
  seurat_obj_dt,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj_dt)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

seurat_obj_dt <- ConstructNetwork(
  seurat_obj_dt,
  tom_name = 'gastruloids',soft_power = 12,corType = "pearson",mergeCutHeight = 0.3,overwrite_tom = TRUE,minModuleSize = 60,deepSplit = 4)

PlotDendrogram(seurat_obj_dt, main='INH hdWGCNA Dendrogram')

seurat_obj_dt <- ModuleEigengenes(
  seurat_obj_dt,
  group.by.vars="orig.ident"
)

# harmonized module eigengenes:
hMEs <- GetMEs(seurat_obj_dt)

# module eigengenes:
MEs <- GetMEs(seurat_obj_dt, harmonized=FALSE)

seurat_obj_dt <- ModuleConnectivity(
  seurat_obj_dt,
  group.by = 'seurat_clusters', group_name = c("0","1","2","3","4","5","6","7","8","9","10","11","12","13")
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj_dt,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

#stitch together with patchwork

wrap_plots(plot_list, ncol=6)

modules <- GetModules(seurat_obj_dt) %>% subset(module == c('pink'))