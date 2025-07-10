# Subset integration on perturbed dataset (iteratively)
ablated3p.subset <- subset_integrate(int_obj = ablated3p.subset, 
                                     perserved_obj = pbmc_3p.ablated, 
                                     subset_ident = "data_5p", 
                                     overlap_clusters = overlapcluster.3p)
ablated5p.subset <- subset_integrate(int_obj = ablated5p.subset, 
                                     perserved_obj = pbmc_5p.ablated, 
                                     subset_ident = "data_3p", 
                                     overlap_clusters = overlapcluster.5p)
# Cluster  + visualization
seurat_visualize_clusters(seurat_obj = ablated3p.subset, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_3p",
                          title = paste0(celltype.keep, " ablated 3p subset integration"))
seurat_visualize_clusters(seurat_obj = ablated5p.subset, 
                          highlight = celltype.keep, 
                          preserve_ident = "data_5p",
                          title = paste0(celltype.keep, " ablated 5p subset integration"))
