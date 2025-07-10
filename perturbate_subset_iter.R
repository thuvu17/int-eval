# Subset integration on perturbed dataset (iteratively)
iter <- iter + 1
if (!is.null(overlapcluster_iter.3p)) {
  ablated3p.subset <- subset_integrate_iter(int_obj = ablated3p.subset, 
                                            perserved_obj = pbmc_3p.ablated, 
                                            subset_ident = "data_5p", 
                                            overlap_clusters = overlapcluster_iter.3p)
  seurat_visualize_clusters(seurat_obj = ablated3p.subset, 
                            highlight = celltype.keep, 
                            preserve_ident = "data_3p",
                            title = paste0(celltype.keep, 
                                           " ablated 3p subset integration ",
                                           "iteration ", iter))
}

if (!is.null(overlapcluster_iter.5p)) {
  ablated5p.subset <- subset_integrate_iter(int_obj = ablated5p.subset, 
                                            perserved_obj = pbmc_5p.ablated, 
                                            subset_ident = "data_3p", 
                                            overlap_clusters = overlapcluster_iter.5p)
  seurat_visualize_clusters(seurat_obj = ablated5p.subset, 
                            highlight = celltype.keep, 
                            preserve_ident = "data_5p",
                            title = paste0(celltype.keep, 
                                           " ablated 5p subset integration ",
                                           "iteration ", iter))
}
