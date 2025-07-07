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
int.balanced <- seurat_clustering(int.balanced)
int.balanced <- seurat_visualize_clusters(int.balanced, 
                                title = "Integrated balanced two batch dataset")

# Marker gene stability analysis
for (celltype.keep in celltype.keep.list) {
  marker_gene_stability(int.balanced, celltype = celltype.keep, title = "balanced",
                        save_path = "results/marker_gene_stability/balanced")
}

# Extract labels and clusters for ARI
int.balanced.label <- int.balanced$anno
int.balanced.cluster <- int.balanced$seurat_clusters
filename <- "balanced"
ari_prep(int.balanced.label, int.balanced.cluster, celltype.keep, filename)
remove(int.balanced.label, int.balanced.cluster, filename)