# Perform perturbation by ablation of all cell types except one

celltype.keep.list <- list("B cell", "NK cell", "Monocyte_FCGR3A", 
                           "Monocyte_CD14", "CD4 T cell", "CD8 T cell")

canonical_markers_list <- list(
  "B cell" = c("MS4A1"),
  "NK cell" = c("GNLY", "NKG7"),
  "Monocyte_FCGR3A" = c("FCGR3A", "MS4A7"),
  "Monocyte_CD14" = c("CD14", "LYZ"),
  "CD4 T cell" = c("IL7R", "CD4", "CCR7"),
  "CD8 T cell" = c("CD8A", "CD8B")
)

# Remove all other cell types
celltype.keep <- celltype.keep.list[[3]]
pbmc_3p.ablated <- subset(pbmc_3p.balanced, subset = celltype == celltype.keep)
pbmc_5p.ablated <- subset(pbmc_5p.balanced, subset = celltype == celltype.keep)

# ==============================================================================
# Perform normal integration
source("~/int-eval/perturbate_v4.R")

# Cluster analysis to identify cluster corresponding to celltype.keep
canonical_markers <- canonical_markers_list[[celltype.keep]]
cluster_analysis(ablated3p, ablated5p, canonical_markers)

# Assign predicted labels by cluster + DEG analysis
pred_clusters.3p <- c(1)
pred_clusters.5p <- c(0)
pred_cells.3p <- WhichCells(ablated3p, expression = seurat_clusters %in% pred_clusters.3p)
pred_cells.5p <- WhichCells(ablated5p, expression = seurat_clusters %in% pred_clusters.5p)
ablated3p$v4_deg_pred <- paste0("non ", celltype.keep)
ablated5p$v4_deg_pred <- paste0("non ", celltype.keep)
ablated3p@meta.data[pred_cells.3p, "v4_deg_pred"] <- celltype.keep
ablated5p@meta.data[pred_cells.5p, "v4_deg_pred"] <- celltype.keep

confusion_matrix(ablated3p, pred_meta = "v4_deg_pred", celltype = celltype.keep,
                 dataset = "Ablated 3p normal Seurat v4")
confusion_matrix(ablated5p, pred_meta = "v4_deg_pred", celltype = celltype.keep,
                 dataset = "Ablated 5p normal Seurat v4")

# Assign predicted labels by overlapping with true labels
pred_clusters.3p <- c(1, 5)
pred_clusters.5p <- c(0)
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

# Marker gene stability
marker_cluster.3p <- c("1", "5")
marker_cluster.5p <- c("0", "1")
marker_gene_stability(seurat_obj = ablated3p, ident1 = marker_cluster.3p, 
                      title = "normal_v4_3p", ident = "seurat_clusters",
                      celltype = celltype.keep,
                      save_path = "results/marker_gene_stability/normal_v4")
marker_gene_stability(seurat_obj = ablated5p, ident1 = marker_cluster.5p, 
                      title = "normal_v4_5p", ident = "seurat_clusters",
                      celltype = celltype.keep,
                      save_path = "results/marker_gene_stability/normal_v4")

# Extract labels and clusters for ARI (NOT FIXED)
# cell_selected <- WhichCells(int.balanced, expression = anno == celltype.keep)
# int.balanced.label <- int.balanced$anno[cell_selected]
# int.balanced.cluster <- int.balanced$seurat_clusters[cell_selected]
# filename <- "balanced"
# ari_prep(int.balanced.label, int.balanced.cluster, celltype.keep, filename)

# ==============================================================================
# Subset integration
overlapcluster.3p <- list(1, 5)
overlapcluster.5p <-list(0, 2, 3, 4)
source("~/int-eval/perturbate_subset.R")

# Marker gene stability
cluster_analysis(ablated3p.subset, ablated5p.subset, canonical_markers,
                 method = "subset integration")
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

# Save results
Idents(ablated3p) <- "pred_anno"
Idents(ablated5p) <- "pred_anno"
marker_gene_stability(seurat_obj = ablated3p, ident = "pred_anno",
                      ident1 = celltype.keep, celltype = celltype.keep, 
                      title = "subset_3p", 
                      save_path = "results/marker_gene_stability/subset")
marker_gene_stability(seurat_obj = ablated5p, ident = "pred_anno",
                      ident1 = celltype.keep, celltype = celltype.keep, 
                      title = "subset_5p", 
                      save_path = "results/marker_gene_stability/subset")
