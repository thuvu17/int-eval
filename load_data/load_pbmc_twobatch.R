library(Seurat)


# Load data
pbmc_3p <- Read10X("data/pbmc8k_filtered_gene_3p")
pbmc_5p <- Read10X_h5("data/pbmc10k_filtered_gene_5p.h5")

# Create Seurat Object
pbmc_3p <- CreateSeuratObject(counts = pbmc_3p, project = "pbmc_3p", min.cells = 3, min.features = 200)
pbmc_5p <- CreateSeuratObject(counts = pbmc_5p[["Gene Expression"]], project = "pbmc_5p", min.cells = 3, min.features = 200)

# Add metadata
pbmc_3p$batch <- "3p"
pbmc_5p$batch <- "5p"

# Add batch name to barcode
pbmc_3p <- RenameCells(pbmc_3p, new.names = paste0("pbmc_3p_", colnames(pbmc_3p)))
pbmc_5p <- RenameCells(pbmc_5p, new.names = paste0("pbmc_5p_", colnames(pbmc_5p)))