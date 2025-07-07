# Perform perturbation by ablation of all cell types except one

celltype.keep.list <- list("Monocyte_CD14", "NK cell", "Monocyte_FCGR3A", 
                           "B cell", "CD4 T cell", "CD8 T cell")

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
marker_gene_stability(seurat_obj = int.balanced, celltype = celltype.keep, 
                      title = "normal_v4",
                      save_path = "results/marker_gene_stability/normal_v4")

# Extract labels and clusters for ARI (NOT FIXED)
cell_selected <- WhichCells(int.balanced, expression = anno == celltype.keep)
int.balanced.label <- int.balanced$anno[cell_selected]
int.balanced.cluster <- int.balanced$seurat_clusters[cell_selected]
filename <- "balanced"
ari_prep(int.balanced.label, int.balanced.cluster, celltype.keep, filename)

# ==============================================================================
# Subset integration
overlapcluster.3p <- list(0, 2, 6)
overlapcluster.5p <-list(0, 1, 6)

ablated3p.subset <- subset_integrate(int_obj = ablated3p, 
                                     perserved_obj = pbmc_3p.ablated, 
                                     subset_ident = "data_5p", 
                                     overlap_clusters = overlapcluster.3p)
ablated5p.subset <- subset_integrate(int_obj = ablated5p, 
                                     perserved_obj = pbmc_5p.ablated, 
                                     subset_ident = "data_3p", 
                                     overlap_clusters = overlapcluster.5p)
# Cluster  + visualization
ablated3p.subset <- seurat_clustering(seurat_obj = abalted3p.subset, resolution = 1.2)
ablated5p.subset <- seurat_clustering(seurat_obj = abalted5p.subset, resolution = 1.2)
seurat_visualize_clusters(seurat_obj = ablated3p.subset, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_3p",
                title = paste0(celltype.keep, " ablated 3p subset integration"))
seurat_visualize_clusters(seurat_obj = ablated5p.subset, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_5p",
                title = paste0(celltype.keep, " ablated 5p subset integration"))

# Marker gene stability (CANNOT USE FINDMARKERS)
marker_gene_stability(seurat_obj = ablated3p.subset, celltype = celltype.keep, 
                      title = "subset", save_path = "results/marker_gene_stability/subset")
marker_gene_stability(seurat_obj = ablated5p.subset, celltype = celltype.keep, 
                      title = "subset", save_path = "results/marker_gene_stability/subset")
