# Perform perturbation by ablation of all cell types except one

celltype.keep.list <- list("B cell", "NK cell", "Monocyte_FCGR3A", 
                           "Monocyte_CD14", "CD4 T cell", "CD8 T cell")

# Remove all other cell types
celltype.keep <- celltype.keep.list[1]
pbmc_3p.ablated <- subset(pbmc_3p.balanced, subset = celltype == celltype.keep)
pbmc_5p.ablated <- subset(pbmc_5p.balanced, subset = celltype == celltype.keep)

# ==============================================================================
# Perform normal integration
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

# Marker gene stability
marker_cluster.3p <- "0"
marker_cluster.5p <- "0"
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
overlapcluster.3p <- list(0, 6)
overlapcluster.5p <-list(0, 7)

ablated3p.subset <- subset_integrate(int_obj = ablated3p, 
                                     perserved_obj = pbmc_3p.ablated, 
                                     subset_ident = "data_5p", 
                                     overlap_clusters = overlapcluster.3p)
ablated5p.subset <- subset_integrate(int_obj = ablated5p, 
                                     perserved_obj = pbmc_5p.ablated, 
                                     subset_ident = "data_3p", 
                                     overlap_clusters = overlapcluster.5p)
# Cluster  + visualization
ablated3p.subset <- seurat_clustering(seurat_obj = ablated3p.subset, resolution = 1.2)
ablated5p.subset <- seurat_clustering(seurat_obj = ablated5p.subset, resolution = 1.2)
seurat_visualize_clusters(seurat_obj = ablated3p.subset, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_3p",
                title = paste0(celltype.keep, " ablated 3p subset integration"))
seurat_visualize_clusters(seurat_obj = ablated5p.subset, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_5p",
                title = paste0(celltype.keep, " ablated 5p subset integration"))
# Get subset labels
ablated3p.chosencells <- WhichCells(ablated3p.subset, expression = Sample != "data_3p")
ablated5p.chosencells <- WhichCells(ablated5p.subset, expression = Sample != "data_5p")
ablated3p$pred_anno <- paste0("non ", celltype.keep)
ablated5p$pred_anno <- paste0("non ", celltype.keep)
ablated3p@meta.data[ablated3p.chosencells, "pred_anno"] <- celltype.keep
ablated5p@meta.data[ablated5p.chosencells, "pred_anno"] <- celltype.keep

# Marker gene stability
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
