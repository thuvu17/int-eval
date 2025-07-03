library(Seurat)
library(patchwork)


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
pbmc_5p.balanced <- downsample(pbmc_5p)

remove(valid_3p, valid_5p, shared_celltypes)

# ==============================================================================
# Unintegrated balanced datasets
unint.balanced <- merge(pbmc_3p.balanced, pbmc_5p.balanced)

# run standard analysis workflow
unint.balanced <- NormalizeData(unint.balanced)
unint.balanced <- FindVariableFeatures(unint.balanced)
unint.balanced <- ScaleData(unint.balanced)
unint.balanced <- RunPCA(unint.balanced)

unint.balanced <- FindNeighbors(unint.balanced, dims = 1:30, reduction = "pca")
unint.balanced <- FindClusters(unint.balanced, resolution = 2, cluster.name = "unintegrated_clusters")

unint.balanced <- RunUMAP(unint.balanced, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(unint.balanced, reduction = "umap.unintegrated", group.by = c("batch", "celltype", "seurat_clusters")) +
  plot_annotation("Unintegrated balanced two batch dataset")

# ==============================================================================
# Integrate the 2 balanced datasets
int.balanced.list <- list(pbmc_3p.balanced, pbmc_5p.balanced)

# normalize and identify variable features for each dataset independently
int.balanced.list <- lapply(X = int.balanced.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# select features that are repeatedly variable across datasets for integration
features <- SelectIntegrationFeatures(object.list = int.balanced.list)

anchors <- FindIntegrationAnchors(object.list = int.balanced.list, 
                                  anchor.features = features)

int.balanced <- IntegrateData(anchorset = anchors)
remove(int.balanced.list, features, anchors)

DefaultAssay(int.balanced) <- "integrated"

# Run the standard workflow for visualization and clustering
int.balanced <- ScaleData(int.balanced, verbose = FALSE)
int.balanced <- RunPCA(int.balanced, npcs = 30, verbose = FALSE)
int.balanced <- RunUMAP(int.balanced, reduction = "pca", dims = 1:30)
int.balanced <- FindNeighbors(int.balanced, reduction = "pca", dims = 1:30)
int.balanced <- FindClusters(int.balanced, resolution = 0.5)

DimPlot(int.balanced, reduction = "umap", group.by = c("batch", "celltype", "seurat_clusters")) +
  plot_annotation("Integrated balanced two batch dataset")
