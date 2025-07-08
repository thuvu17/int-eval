# Subset integration on perturbed dataset (run this after perturbate_v4)
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
