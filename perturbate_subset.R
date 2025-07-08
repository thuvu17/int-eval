# Subset integration on perturbed dataset (run this after perturbate_v4)
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
