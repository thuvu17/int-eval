# Perform perturbation by ablation of all cell types except one

celltype.keep.list <- list("Monocyte_CD14", "NK cell", "Monocyte_FCGR3A", 
                           "B cell", "CD4 T cell", "CD8 T cell")

# Remove all other cell types
celltype.keep <- celltype.keep.list[1]
pbmc_3p.ablated <- subset(pbmc_3p.balanced, subset = celltype.keep)
pbmc_5p.ablated <- subset(pbmc_5p.balanced, subset = celltype.keep)

# ==============================================================================
# Perform normal integration
ablated3p.list <- list(pbmc_3p.ablated, pbmc_5p.balanced)
ablated5p.list <- list(pbmc_5p.ablated, pbmc_3p.balanced)

ablated3p <- seuratv4_integrate(ablated3p.list)
ablated5p <- seuratv4_integrate(ablated5p.list)

remove(ablated3p.list, ablated5p.list)

# ==============================================================================
# Subsetting integration

