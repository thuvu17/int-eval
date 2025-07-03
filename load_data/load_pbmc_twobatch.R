library(zellkonverter)
library(Seurat)
library(SingleCellExperiment)


# Set file path
file_path.1 <- "data/pbmc_2_batch_base_balanced/tran_exp5_pbmc_batch1_balanced.h5ad"
file_path.2 <- "data/pbmc_2_batch_base_balanced/tran_exp5_pbmc_batch2_balanced.h5ad"

# Read h5ad files
pbmc_3p.balanced <- readH5AD(file_path.1)
pbmc_5p.balanced <- readH5AD(file_path.2)

# Create Seurat Object
pbmc_3p.balanced <- as.Seurat(pbmc_3p.balanced, counts = "X", data = "X")
pbmc_5p.balanced <- as.Seurat(pbmc_5p.balanced, counts = "X", data = "X")
