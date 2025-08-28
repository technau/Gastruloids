#This script takes you through projecting the module gene scores from the script hdWGCNA_treatment response.R onto the query gastruloid dataset 
# Load required libraries
library(Seurat)
library(ggplot2)
library(dplyr)
library(ggpubr)

load('workspace from hdWGCNA_treatment response.R')

# ------------------------
# Step 1: Define Module Genes
# ------------------------
module_genes <- modules$gene_name  # replace with your genes
module_genes <- module_genes[module_genes %in% rownames(seurat_obj_dt)]

# ------------------------
# Step 2: Add Module Score to Seurat Object
# ------------------------
load('object from scRNAseq analysis gastruloids.R AFTER SingleR_labeltransfer.R')
seurat_obj_g <- AddModuleScore(
  object = seurat_obj_g,
  features = list(module_genes),
  name = "Module1"
)

# ------------------------
# Step 3: Prepare Metadata for Plotting
# ------------------------
metadata_df <- seurat_obj_g@meta.data %>%
  dplyr::mutate(module_score = Module11) %>%
  dplyr::select(singleR.labels, newIdents, module_score)

# Optional: Filter only annotated cells
metadata_df <- metadata_df %>%
  dplyr::filter(!is.na(singleR.labels), !is.na(newIdents))

# ------------------------
# Step 4: Violin Plot Grouped by Celltype and Condition
# ------------------------
# Combine the two metadata columns into one grouping variable
metadata_df$Group <- paste(metadata_df$singleR.labels, metadata_df$newIdents, sep = "_")
metadata_df$Group <- factor(
  metadata_df$Group,
  levels = c("Endoderm_wt.24h", "Endoderm_agg.4h",
             "Endoderm_agg.12h", "Endoderm_agg.24h",
             "Endoderm_agg.48h","Ectoderm_wt.24h", "Ectoderm_agg.4h",
             "Ectoderm_agg.12h", "Ectoderm_agg.24h",
             "Ectoderm_agg.48h","Mesoderm_wt.24h", "Mesoderm_agg.4h",
             "Mesoderm_agg.12h", "Mesoderm_agg.24h",
             "Mesoderm_agg.48h","Neuronal_wt.24h", "Neuronal_agg.4h",
             "Neuronal_agg.12h", "Neuronal_agg.24h",
             "Neuronal_agg.48h","Cnidocytes_wt.24h", "Cnidocytes_agg.4h",
             "Cnidocytes_agg.12h", "Cnidocytes_agg.24h",
             "Cnidocytes_agg.48h","Secretory Progenitor_wt.24h", "Secretory Progenitor_agg.4h",
             "Secretory Progenitor_agg.12h", "Secretory Progenitor_agg.24h",
             "Secretory Progenitor_agg.48h","Gland Cells_wt.24h", "Gland Cells_agg.4h",
             "Gland Cells_agg.12h", "Gland Cells_agg.24h",
             "Gland Cells_agg.48h"))
# Plot
ggplot(metadata_df, aes(x = Group, y = module_score, fill = newIdents)) +
  geom_violin(scale = "width", trim = FALSE) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  ylab("Module Score") +
  xlab("Cell Type and Timepoints") +
  ggtitle("Pink module Score Across Cell Types and Timepoints") +
  scale_fill_manual(values = c("wt.24h" = "#56B4E9", "agg.4h" = "#D55E00","agg.12h" = "black","agg.24h" = "pink","agg.48h" = "red"))

#Use differentially expressed genes specific only to the mesodermal cluster and treat is as a mesodermal module
metadata_df1$module_name <- "Pink module"
metadata_df2$module_name <- "Mesodermal markers"  # Replace with actual name

# Rename columns if needed
names(metadata_df1)[names(metadata_df1) == "score_column_name"] <- "module_score"
names(metadata_df2)[names(metadata_df2) == "score_column_name"] <- "module_score"

combined_df <- rbind(metadata_df2, metadata_df1)

ggplot(combined_df, aes(x = Group, y = module_score, fill = newIdents)) +
  geom_violin(scale = "width", trim = FALSE) +
  facet_wrap(~ module_name, ncol = 1, scales = "free_y") +  # One column, separate y-axes
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  xlab("Cell Type and Timepoints") +
  ylab("Module Score") +
  scale_fill_manual(values = c("wt.24h" = "#56B4E9", 
                               "agg.4h" = "#D55E00",
                               "agg.12h" = "black",
                               "agg.24h" = "pink",
                               "agg.48h" = "red")) +
  ggtitle("Module Scores Across Cell Types and Timepoints")

#Ranking the module scores

library(dplyr)

# Step 1: Create Group variable (cell type + timepoint)
seurat_obj_g$Group <- paste(seurat_obj_g$singleR.labels, seurat_obj_g$newIdents, sep = "_")

# Step 2: Compute average module score per (cell type × timepoint)
pseudobulk_df <- seurat_obj_g@meta.data %>%
  mutate(module_score = seurat_obj_g$Module11) %>%
  group_by(singleR.labels, newIdents) %>%
  summarise(avg_score = mean(module_score), .groups = "drop")

# Step 3: Rank module scores within each timepoint across cell types
ranked_scores <- pseudobulk_df %>%
  group_by(newIdents) %>%
  arrange(desc(avg_score)) %>%
  mutate(rank = rank(-avg_score, ties.method = "first")) %>%  # Rank in descending order (highest = rank 1)
  ungroup()

# View result
head(ranked_scores)
