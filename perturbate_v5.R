# Perform normal integration on perturbed datasets

ablated3p.v5 <- merge(pbmc_3p.ablated, pbmc_5p.balanced)
ablated5p.v5 <- merge(pbmc_5p.ablated, pbmc_3p.balanced)

ablated3p.harmony <- seuratv5_integrate(ablated3p.v5, method = "harmony")
ablated5p.harmony <- seuratv5_integrate(ablated5p.v5, method = "harmony")

ablated3p.scvi <- seuratv5_integrate(ablated3p.v5, method = "scvi")
ablated5p.scvi <- seuratv5_integrate(ablated5p.v5, method = "scvi")

# Standard workflow for visualization and clustering
ablated3p.harmony <- seurat_clustering(ablated3p.harmony, reduction = "harmony")
ablated5p.harmony <- seurat_clustering(ablated5p.harmony, reduction = "harmony")

ablated3p.scvi <- seurat_clustering(ablated3p.scvi, reduction = "integrated.scvi")
ablated5p.scvi <- seurat_clustering(ablated5p.scvi, reduction = "integrated.scvi")

seurat_visualize_clusters(seurat_obj = ablated3p.harmony, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_3p",
                          title = paste0(celltype.keep, " ablated 3p Harmony"))
seurat_visualize_clusters(seurat_obj = ablated5p.harmony, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_5p",
                          title = paste0(celltype.keep, " ablated 5p Harmony"))
seurat_visualize_clusters(seurat_obj = ablated3p.scvi, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_3p",
                          title = paste0(celltype.keep, " ablated 3p scVI"))
seurat_visualize_clusters(seurat_obj = ablated5p.scvi, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_5p",
                          title = paste0(celltype.keep, " ablated 5p scVI"))

# ==============================================================================
# Cluster analysis to identify cluster corresponding to celltype.keep
canonical_markers <- canonical_markers_list[[celltype.keep]]
ablated3p.harmony <- JoinLayers(ablated3p.harmony)
ablated5p.harmony <- JoinLayers(ablated5p.harmony)
deg_results <- cluster_analysis(ablated3p.harmony, ablated5p.harmony, 
                                markers = canonical_markers, top_n = 50)

# ------------------------------------------------------------------------------
# Assign predicted labels by cluster + DEG analysis
pred_clusters.3p <- deg_results$clusters.3p
pred_clusters.5p <- deg_results$clusters.5p
pred_cells.3p <- WhichCells(ablated3p.harmony, expression = seurat_clusters %in% pred_clusters.3p)
pred_cells.5p <- WhichCells(ablated5p.harmony, expression = seurat_clusters %in% pred_clusters.5p)
ablated3p.harmony$harmony_deg_pred <- paste0("non ", celltype.keep)
ablated5p.harmony$harmony_deg_pred <- paste0("non ", celltype.keep)
ablated3p.harmony@meta.data[pred_cells.3p, "harmony_deg_pred"] <- celltype.keep
ablated5p.harmony@meta.data[pred_cells.5p, "harmony_deg_pred"] <- celltype.keep

confusion_matrix(ablated3p.harmony, pred_meta = "harmony_deg_pred", celltype = celltype.keep,
                 dataset = "Ablated 3p Harmony v5 + deg")
confusion_matrix(ablated5p.harmony, pred_meta = "harmony_deg_pred", celltype = celltype.keep,
                 dataset = "Ablated 5p Harmony v5 + deg")

# ------------------------------------------------------------------------------
# Assign predicted labels by overlapping with true labels
pred_clusters.3p <- c(0, 10)
pred_clusters.5p <- c(0, 8)
pred_cells.3p <- WhichCells(ablated3p, expression = seurat_clusters %in% pred_clusters.3p)
pred_cells.5p <- WhichCells(ablated5p, expression = seurat_clusters %in% pred_clusters.5p)
pred_cells.3p <- WhichCells(ablated3p.harmony, expression = seurat_clusters %in% pred_clusters.3p)
pred_cells.5p <- WhichCells(ablated5p.harmony, expression = seurat_clusters %in% pred_clusters.5p)
ablated3p.harmony$harmony_deg_pred <- paste0("non ", celltype.keep)
ablated5p.harmony$harmony_deg_pred <- paste0("non ", celltype.keep)
ablated3p.harmony@meta.data[pred_cells.3p, "harmony_overlap_pred"] <- celltype.keep
ablated5p.harmony@meta.data[pred_cells.5p, "harmony_overlap_pred"] <- celltype.keep

confusion_matrix(ablated3p.harmony, pred_meta = "harmony_overlap_pred", celltype = celltype.keep,
                 dataset = "Ablated 3p Harmony v5 + overlap")
confusion_matrix(ablated5p.harmony, pred_meta = "harmony_overlap_pred", celltype = celltype.keep,
                 dataset = "Ablated 5p Harmony v5 + overlap")
