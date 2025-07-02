# To create balanced two batch PBMC, downsample every cell type cluster to 200 cells
# First, only keep cell types with >= 200 cells
valid_3p <- names(which(table(pbmc_3p$celltype) >= 200))
valid_5p <- names(which(table(pbmc_5p$celltype) >= 200))
shared_celltypes <- intersect(valid_3p, valid_5p)

pbmc_3p <- subset(pbmc_3p, subset = celltype %in% shared_celltypes)
pbmc_5p <- subset(pbmc_5p, subset = celltype %in% shared_celltypes)