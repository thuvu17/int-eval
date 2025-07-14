# Perform perturbation by ablation of all cell types except one

celltype.keep.list <- list("B cell", "NK cell", "Monocyte_FCGR3A", 
                           "Monocyte_CD14", "CD4 T cell", "CD8 T cell")

canonical_markers_list <- list(
  "B cell" = c("MS4A1"),
  "NK cell" = c(" ", "NKG7"),
  "Monocyte_FCGR3A" = c("FCGR3A", "MS4A7"),
  "Monocyte_CD14" = c("CD14", "LYZ"),
  "CD4 T cell" = c("IL7R", "CD4", "CCR7"),
  "CD8 T cell" = c("CD8A", "CD8B")
)

# Remove all other cell types
celltype.keep <- celltype.keep.list[[2]]
pbmc_3p.ablated <- subset(pbmc_3p.balanced, subset = celltype == celltype.keep)
pbmc_5p.ablated <- subset(pbmc_5p.balanced, subset = celltype == celltype.keep)

# ==============================================================================
# Perform normal integration
source("perturbate_v4.R")

# Cluster analysis to identify cluster corresponding to celltype.keep
canonical_markers <- canonical_markers_list[[celltype.keep]]
cluster_analysis(ablated3p, ablated5p, canonical_markers)

# Assign predicted labels by cluster + DEG analysis
pred_clusters.3p <- c(0, 5)
pred_clusters.5p <- c(4, 5)
pred_cells.3p <- WhichCells(ablated3p, expression = seurat_clusters %in% pred_clusters.3p)
pred_cells.5p <- WhichCells(ablated5p, expression = seurat_clusters %in% pred_clusters.5p)
ablated3p$v4_deg_pred <- paste0("non ", celltype.keep)
ablated5p$v4_deg_pred <- paste0("non ", celltype.keep)
ablated3p@meta.data[pred_cells.3p, "v4_deg_pred"] <- celltype.keep
ablated5p@meta.data[pred_cells.5p, "v4_deg_pred"] <- celltype.keep

confusion_matrix(ablated3p, pred_meta = "v4_deg_pred", celltype = celltype.keep,
                 dataset = "Ablated 3p normal Seurat v4 + deg")
confusion_matrix(ablated5p, pred_meta = "v4_deg_pred", celltype = celltype.keep,
                 dataset = "Ablated 5p normal Seurat v4 + deg")

# Assign predicted labels by overlapping with true labels
pred_clusters.3p <- c(0, 9)
pred_clusters.5p <- c(0, 8)
pred_cells.3p <- WhichCells(ablated3p, expression = seurat_clusters %in% pred_clusters.3p)
pred_cells.5p <- WhichCells(ablated5p, expression = seurat_clusters %in% pred_clusters.5p)
ablated3p$v4_overlap_pred <- paste0("non ", celltype.keep)
ablated5p$v4_overlap_pred <- paste0("non ", celltype.keep)
ablated3p@meta.data[pred_cells.3p, "v4_overlap_pred"] <- celltype.keep
ablated5p@meta.data[pred_cells.5p, "v4_overlap_pred"] <- celltype.keep

confusion_matrix(ablated3p, pred_meta = "v4_overlap_pred", celltype = celltype.keep,
                 dataset = "Ablated 3p normal Seurat v4 + overlap")
confusion_matrix(ablated5p, pred_meta = "v4_overlap_pred", celltype = celltype.keep,
                 dataset = "Ablated 5p normal Seurat v4 + overlap")

# ==============================================================================
# Subset integration
overlapcluster.3p <- list(0, 9)
overlapcluster.5p <-list(0, 8)
source("perturbate_subset.R")

# Subset iteratively
iter <- 1
overlapcluster_iter.3p <- NULL
overlapcluster_iter.5p <- list(0, 2, 3, 4, 5)
source("perturbate_subset_iter.R")
length(WhichCells(ablated5p.subset, expression = celltype == "Monocyte_FCGR3A"))

# Get subset labels
pred_cells.3p <- WhichCells(ablated3p.subset)
pred_cells.5p <- WhichCells(ablated5p.subset)
ablated3p$subset_pred <- paste0("non ", celltype.keep)
ablated5p$subset_pred <- paste0("non ", celltype.keep)
ablated3p@meta.data[pred_cells.3p, "subset_pred"] <- celltype.keep
ablated5p@meta.data[pred_cells.5p, "subset_pred"] <- celltype.keep

confusion_matrix(ablated3p, pred_meta = "subset_pred", celltype = celltype.keep,
                 dataset = "Ablated 3p subset integration")
confusion_matrix(ablated5p, pred_meta = "subset_pred", celltype = celltype.keep,
                 dataset = "Ablated 5p subset integration")

