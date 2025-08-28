#Script to generate reference atlas object for annotating query cell clusters 
#you can download the Nematostella NV2 atlas object from https://cells.ucsc.edu/?ds=sea-anemone-atlas+Nv2+all
load("AllData.Robj")
Idents(AllData) <- "orig.ident"

#subsetting only relevant early embryonic stages to annotate early embryonic gastruloids
reference <- subset(AllData,idents = c("D.2d planula","D.2d2 planula","D.3d planula","D.4d planula"))
wt.24h[["RNA"]] <- JoinLayers(wt.24h[["RNA"]])
reference <- merge(x = wt.24h,y = list(reference))
reference[["RNA"]] <- JoinLayers(reference[["RNA"]])
assay5 <- as(reference[["RNA"]], Class = "Assay5")

reference5 <- CreateSeuratObject(assay5, meta.data = reference@meta.data)
reference5 <- NormalizeData(reference5)
reference5 <- FindVariableFeatures(reference5)
reference5 <- ScaleData(reference5)
reference5 <- RunPCA(reference5)
reference5[["RNA"]] <- split(reference5[["RNA"]], f = reference5$orig.ident)
reference5 <- IntegrateLayers(object = reference5, method = HarmonyIntegration, orig.reduction = "pca",new.reduction = "harmony")
reference5[["RNA"]] <- JoinLayers(reference5[["RNA"]])
reference5 <- FindNeighbors(reference5, reduction = "harmony", dims = 1:30)
reference5 <- FindClusters(reference5, resolution = 0.4)
reference5 <- RunUMAP(reference5,dims = 1:30,reduction = "harmony",return.model = TRUE,umap.method = "uwot-learn")
new.cluster.ids <- c(`0` = "Ectoderm",`1` = "Ectoderm",`2` = "Ectoderm",`3` = "Ectoderm",`4` = "Ectoderm",`5` = "Secretory Progenitor",`6` = "Ectoderm",`7` = "Mesoderm",`8` = "Mesoderm",`9` = "Cnidocytes",`10` = "Endoderm",`11` = "Endoderm",`12` = "Gland Cells",`13` = "Neuronal",`14` = "Neuronal",`15` = "Cnidocytes",`16` = "Cnidocytes",`17` = "Gland Cells",`18` = "Neuronal")
c7 <- c(`Ectoderm` = "dodgerblue2", `Endoderm` = "green4",`Mesoderm` = "#E31A1C",`Secretory Progenitor` = "#FF7F00",`Cnidocytes` = "black", `Gland Cells` = "gold1",`Neuronal` = "orchid1")
c20_native <- c(`0` = "dodgerblue2", `1` = "#E31A1C",`2` = "green4",`3` = "#6A3D9A",`4` = "#FF7F00",`5` = "black", `6` = "gold1",`7` = "skyblue2",`8` = "palegreen2",`9` = "gray70", `10` = "khaki2",`11` = "maroon", `12` = "orchid1", `13` = "deeppink1",`14` =  "blue1", `15` = "steelblue4",`16` = "darkturquoise",`17` =  "green1", `18` = "yellow4", `19` = "brown",`20` = "yellow",`21` = "turquoise",`22` = "orange")
names(new.cluster.ids) <- levels(reference5)
reference5 <- RenameIdents(reference5, new.cluster.ids)
reference5$newIDs <- reference5@active.ident
Idents(reference5) <- "newIDs"

save(file = "Reference.Robj")