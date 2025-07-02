# To create balanced two batch PBMC, downsample every cell type cluster to 200 cells
# First, only keep cell types with >= 200 cells
valid_3p <- names(which(table(pbmc_3p$celltype) >= 200))
valid_5p <- names(which(table(pbmc_5p$celltype) >= 200))
shared_celltypes <- intersect(valid_3p, valid_5p)

pbmc_3p <- subset(pbmc_3p, subset = celltype %in% shared_celltypes)
pbmc_5p <- subset(pbmc_5p, subset = celltype %in% shared_celltypes)

# Then, randomly select only 200 cells in each cell type
downsample <- function(seurat_obj, group_var = "celltype", target_n = 200) {
  celltypes <- table(seurat_obj[[group_var]][, 1])
  valid_types <- names(celltypes[celltypes >= target_n])
  
  selected_cells <- unlist(lapply(valid_types, function(ct) {
    cells <- rownames(seurat_obj@meta.data[seurat_obj@meta.data[[group_var]] == ct, ])
    sample(cells, target_n)
  }))
  
  subset(seurat_obj, cells = selected_cells)
}

pbmc_3p.balanced <- downsample(pbmc_3p)
pbmc_5p_balanced <- downsample(pbmc_5p)
