# Perform normal integration on perturbed datasets

ablated3p.list <- list(pbmc_3p.ablated, pbmc_5p.balanced)
ablated5p.list <- list(pbmc_5p.ablated, pbmc_3p.balanced)

ablated3p <- seuratv4_integrate(ablated3p.list)
ablated5p <- seuratv4_integrate(ablated5p.list)

# Standard workflow for visualization and clustering
ablated3p <- seurat_clustering(ablated3p)
ablated5p <- seurat_clustering(ablated5p)
seurat_visualize_clusters(seurat_obj = ablated3p, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_3p",
                          title = paste0(celltype.keep, " ablated 3p normal integration"))
seurat_visualize_clusters(seurat_obj = ablated5p, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_5p",
                          title = paste0(celltype.keep, " ablated 5p normal integration"))
remove(ablated3p.list, ablated5p.list)

# ==============================================================================
# Cluster analysis to identify cluster corresponding to celltype.keep
canonical_markers <- canonical_markers_list[[celltype.keep]]
cluster_analysis(ablated3p, ablated5p, canonical_markers)

# ------------------------------------------------------------------------------
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

# ------------------------------------------------------------------------------
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