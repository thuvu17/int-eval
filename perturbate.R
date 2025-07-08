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
celltype.keep <- celltype.keep.list[[1]]
pbmc_3p.ablated <- subset(pbmc_3p.balanced, subset = celltype == celltype.keep)
pbmc_5p.ablated <- subset(pbmc_5p.balanced, subset = celltype == celltype.keep)

# ==============================================================================
# Perform normal integration
source("~/int-eval/perturbate_v4.R")

# Cluster analysis to identify cluster corresponding to celltype.keep
canonical_markers <- canonical_markers_list[[celltype.keep]]
cluster_analysis(ablated3p, ablated5p, canonical_markers)

# Marker gene stability
marker_cluster.3p <- c("1", "5")
marker_cluster.5p <- c("0", "6")
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
overlapcluster.5p <-list(0, 6)
source("~/int-eval/perturbate_subset.R")

# Marker gene stability
cluster_analysis(ablated3p.subset, ablated5p.subset, canonical_markers)
# Get subset labels
ablated3p.chosencells <- WhichCells(ablated3p.subset, expression = Sample != "data_3p")
ablated5p.chosencells <- WhichCells(ablated5p.subset, expression = Sample != "data_5p")
ablated3p$pred_anno <- paste0("non ", celltype.keep)
ablated5p$pred_anno <- paste0("non ", celltype.keep)
ablated3p@meta.data[ablated3p.chosencells, "pred_anno"] <- celltype.keep
ablated5p@meta.data[ablated5p.chosencells, "pred_anno"] <- celltype.keep

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
