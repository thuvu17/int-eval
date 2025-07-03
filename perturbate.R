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
ablated3p <- seurat_visualize_clusters(ablated3p, highlight = celltype.keep,
               title = paste0(celltype.keep, " ablated 3p normal integration"))
ablated5p <- seurat_visualize_clusters(ablated5p, highlight = celltype.keep,
               title = paste0(celltype.keep, " ablated 5p normal integration"))

remove(ablated3p.list, ablated5p.list)

# ==============================================================================
# Subset integration

