library(Seurat)
library(patchwork)


# Unintegrated balanced datasets
unint.balanced <- merge(pbmc_3p.balanced, pbmc_5p.balanced)

# run standard analysis workflow
unint.balanced <- NormalizeData(unint.balanced)
unint.balanced <- FindVariableFeatures(unint.balanced)
unint.balanced <- ScaleData(unint.balanced)
unint.balanced <- RunPCA(unint.balanced)

unint.balanced <- FindNeighbors(unint.balanced, dims = 1:30, reduction = "pca")
unint.balanced <- FindClusters(unint.balanced, resolution = 2, 
                               cluster.name = "unintegrated_clusters")

unint.balanced <- RunUMAP(unint.balanced, dims = 1:30, reduction = "pca", 
                          reduction.name = "umap.unintegrated")
DimPlot(unint.balanced, reduction = "umap.unintegrated", 
        group.by = c("batch", "celltype", "seurat_clusters")) +
  plot_annotation("Unintegrated balanced two batch dataset")

# ==============================================================================
# Integrate the 2 balanced datasets
int.balanced.list <- list(pbmc_3p.balanced, pbmc_5p.balanced)
int.balanced <- seuratv4_integrate(int.balanced.list)
remove(int.balanced.list)

# Run the standard workflow for visualization and clustering
int.balanced <- ScaleData(int.balanced, verbose = FALSE)
int.balanced <- RunPCA(int.balanced, npcs = 30, verbose = FALSE)
int.balanced <- RunUMAP(int.balanced, reduction = "pca", dims = 1:30)
int.balanced <- FindNeighbors(int.balanced, reduction = "pca", dims = 1:30)
int.balanced <- FindClusters(int.balanced, resolution = 0.5)

DimPlot(int.balanced, reduction = "umap", 
        group.by = c("batch", "celltype", "seurat_clusters")) +
  plot_annotation("Integrated balanced two batch dataset")
