library(Seurat)
library(ggplot2)
library(patchwork)


# Perform Seurat v4 integration
seuratv4_integrate <- function(obj.list, nfeatures = 2000, anchor.features = NULL) {
  # Normalize and find variable features
  obj.list <- lapply(obj.list, function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = nfeatures)
    return(x)
  })
  
  # Find anchors and integrate
  features <- SelectIntegrationFeatures(object.list = obj.list)
  anchors <- FindIntegrationAnchors(object.list = obj.list, 
                                    anchor.features = features)
  integrated <- IntegrateData(anchorset = anchors)
  DefaultAssay(integrated) <- "integrated"
  
  return(integrated)
}

# Run standard workflow for visualization and clustering
seurat_visualize_clusters <- function(seurat_obj, dims = 1:30, highlight = NULL,
                                      resolution = 0.5, title = "") {
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = max(dims), verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = "pca", dims = dims)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)
  
  p <- DimPlot(seurat_obj, reduction = "umap", 
               group.by = c("batch", "celltype", "seurat_clusters"))
  if (!is.null(highlight)) {
    highlighted_cells <- WhichCells(seurat_obj, expression = celltype == highlight)
    p.highlight <- DimPlot(seurat_obj, reduction = "umap", cells.highlight = highlighted_cells) +
      ggtitle(paste("Highlight:", highlight))
    print(p + p.highlight + plot_annotation(title))
  } else {
    print(p + plot_annotation(title))
  }
  
  return(seurat_obj)
}

# Top 10 marker genes stability analysis
marker_gene_stability <- function(seurat_obj, celltype, title, save_path = NULL) {
  Idents(seurat_obj) <- "celltype"
  markers <- FindConservedMarkers(seurat_obj, ident.1 = celltype, 
                                  grouping.var = "batch", assay = "originalexp",
                                  verbose = FALSE)
  
  cat(paste0(celltype, " top 10 marker genes (", title, "):\n"))
  print(head(markers, 10))
  
  if (!is.null(save_path)) {
    filename <- paste0(celltype, "_", title)
    saveRDS(markers, file = paste0(save_path, "/", filename, ".rds"))
    write.csv(markers, file = paste0(save_path, "/", filename, ".csv"))
  }
}
