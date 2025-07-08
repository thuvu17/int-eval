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

