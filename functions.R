library(Seurat)
library(SeuratWrappers)
library(ggplot2)
library(patchwork)
library(dplyr)


# Perform Seurat v4 integration
seuratv4_integrate <- function(obj.list, nfeatures = 2000) {
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


# Perform Seurat v5 integration (Harmony or scVI)
seuratv5_integrate <- function(obj, method) {
  # Normalize and find variable features
  obj[["originalexp"]] <- split(obj[["originalexp"]], f = obj$Sample)
  obj <- NormalizeData(obj)
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj)
  
  obj <- FindNeighbors(obj, dims = 1:30, reduction = "pca")
  obj <- FindClusters(obj, resolution = 2, cluster.name = "unintegrated_clusters")
  obj <- RunUMAP(obj, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
  
  if (method == "harmony") {
    obj <- IntegrateLayers(
      object = obj, method = HarmonyIntegration,
      orig.reduction = "pca", new.reduction = "harmony",
      verbose = FALSE
    )
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30)
    obj <- FindClusters(obj, resolution = 2, cluster.name = "harmony_clusters")
    obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
  } 
  else if (method == "scvi") {
      obj <- IntegrateLayers(
        object = obj, method = scVIIntegration,
        new.reduction = "integrated.scvi",
        conda_env = "C:/Users/vnmt1/miniconda3/envs/scvi-env", verbose = FALSE
      )
      obj <- FindNeighbors(obj, reduction = "integrated.scvi", dims = 1:30)
      obj <- FindClusters(obj, resolution = 2, cluster.name = "scvi_clusters")
      obj <- RunUMAP(obj, reduction = "integrated.scvi", dims = 1:30, reduction.name = "umap.cca")
  }
  
  return(obj)
}


# Run standard workflow for visualization and clustering
seurat_clustering <- function(seurat_obj, dims = 1:30, resolution = 1.0, reduction = "pca") {
  seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, npcs = max(dims), verbose = FALSE)
  seurat_obj <- RunUMAP(seurat_obj, reduction = reduction, dims = dims)
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction, dims = dims)
  seurat_obj <- FindClusters(seurat_obj, resolution = resolution)

  return(seurat_obj)
}


# Visualize clusters
seurat_visualize_clusters <- function(seurat_obj, highlight = NULL, 
                                      title = "", preserve_ident = NULL) {
  p <- DimPlot(seurat_obj, reduction = "umap", 
               group.by = c("batch", "celltype", "seurat_clusters"))
  if (!is.null(highlight)) {
    highlighted_cells <- WhichCells(seurat_obj, 
                                    expression = Sample == preserve_ident)
    p.highlight <- DimPlot(seurat_obj, reduction = "umap", 
                           cells.highlight = highlighted_cells) +
      ggtitle("Highlight ablated")
    print(p + p.highlight + plot_annotation(title))
  } else {
    print(p + plot_annotation(title))
  }
}


# DEG analysis: Get top 20 marker genes by p-val and log2fc
deg_analysis <- function(seurat_obj, n = 20) {
  DefaultAssay(seurat_obj) <- "originalexp"
  markers <- FindAllMarkers(seurat_obj, only.pos = TRUE)
  top20.markers <- markers %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    dplyr::filter(p_val_adj < 0.05) %>%
    slice_min(order_by = p_val_adj, n = n, with_ties = FALSE)
  
  return(top20.markers)
}


# Cluster analysis
cluster_analysis <- function(seurat_obj.3p, seurat_obj.5p, markers,
                             method = "normal integration", top_n = 20) {
  # DefaultAssay(seurat_obj.3p) <- "originalexp"
  # DefaultAssay(seurat_obj.5p) <- "originalexp"
  # Print plots for 3p
  # featureplot.3p <- FeaturePlot(seurat_obj.3p, features = markers)
  # clusterplot.3p <- DimPlot(seurat_obj.3p, reduction = "umap", group.by = c("seurat_clusters"))
  # title.3p <- paste0(celltype.keep, " ablated 3p ", method)
  # print(featureplot.3p + clusterplot.3p + plot_annotation(title.3p))
  # Print plots for 5p
  # featureplot.5p <- FeaturePlot(seurat_obj.5p, features = markers)
  # clusterplot.5p <- DimPlot(seurat_obj.5p, reduction = "umap", group.by = c("seurat_clusters"))
  # title.5p <- paste0(celltype.keep, " ablated 5p ", method)
  # print(featureplot.5p + clusterplot.5p + plot_annotation(title.5p))
  # DEG analysis
  results.3p <- deg_analysis(seurat_obj.3p, top_n)
  results.5p <- deg_analysis(seurat_obj.5p, top_n)
  
  valid_clusters_3p <- c()
  valid_clusters_5p <- c()
  
  # Print top DEGs for 3p
  cat("\n===== Top 20 DEGs for", celltype.keep, "in 3p =====\n")
  for (clust in unique(results.3p$cluster)) {
    cat("\nCluster", clust, ":\n")
    top_genes <- results.3p %>% dplyr::filter(cluster == clust)
    print(top_genes)
    if (all(markers %in% top_genes$gene)) {
      valid_clusters_3p <- c(valid_clusters_3p, clust)
      print("Cluster chosen")
    }
  }
  # Print top DEGs for 5p
  cat("\n===== Top 20 DEGs for", celltype.keep, "in 5p =====\n")
  for (clust in unique(results.5p$cluster)) {
    cat("\nCluster", clust, ":\n")
    top_genes <- results.5p %>% dplyr::filter(cluster == clust)
    print(top_genes)
    if (all(markers %in% top_genes$gene)) {
      valid_clusters_5p <- c(valid_clusters_5p, clust)
      print("Cluster chosen")
    }
  }
  
  return(list(clusters.3p = valid_clusters_3p, clusters.5p = valid_clusters_5p))
}


# Top 10 marker genes stability analysis
marker_gene_stability <- function(seurat_obj, ident = "celltype", celltype, 
                                  ident1, title, save_path = NULL) {
  Idents(seurat_obj) <- ident
  markers <- FindMarkers(seurat_obj, ident.1 = ident1, assay = "originalexp",
                         verbose = FALSE)
  
  # Top 10 markers by smallest p values
  top10.minpval <- markers %>%
    filter(avg_log2FC > 1) %>%
    slice_min(order_by = p_val_adj, n = 10, with_ties = FALSE)
  # Top 10 markers by log fold change
  top10.logfc <- markers %>%
    filter(avg_log2FC > 1) %>%
    slice_max(order_by = avg_log2FC, n = 10, with_ties = FALSE)
  
  # Print to console
  cat(paste0("\n", celltype, " top 10 marker genes by adj. p-value (", title, "):\n"))
  print(top10.minpval)
  cat(paste0("\n", celltype, " top 10 marker genes by log2FC (", title, "):\n"))
  print(top10.logfc)
  
  # Save data
  if (!is.null(save_path)) {
    filename <- paste0(celltype, "_", title)
    saveRDS(markers, file = paste0(save_path, "/", filename, ".rds"))
    write.csv(top10.minpval, file = paste0(filename, "_top10_minpval.csv"))
    write.csv(top10.logfc, file = paste0(filename, "_top10_logfc.csv"))
  }
}


# Perform subset integration
subset_integrate <- function(int_obj, perserved_obj, subset_ident, overlap_clusters) {
  # Subset overlapping cells
  overlap_cells <- WhichCells(int_obj, expression = Sample == subset_ident & 
                                seurat_clusters %in% overlap_clusters)
  # Extract counts layer and metadata
  counts <- GetAssayData(int_obj, assay = "originalexp", slot = "counts")[, overlap_cells]
  meta <- int_obj@meta.data[overlap_cells, ]
  # Create new Seurat object
  overlap_obj <- CreateSeuratObject(counts, meta.data = meta)
  
  # Perform integration with subset
  subset.list <- c(perserved_obj, overlap_obj)
  subsetint_obj <- seuratv4_integrate(subset.list)
  subsetint_obj <- seurat_clustering(subsetint_obj)
  
  return(subsetint_obj)
}

subset_integrate_iter <- function(int_obj, perserved_obj, subset_ident, overlap_clusters) {
  # Subset overlapping cells
  overlap_cells <- WhichCells(int_obj, expression = Sample == subset_ident & 
                                seurat_clusters %in% overlap_clusters)
  # Extract counts layer and metadata
  counts <- GetAssayData(int_obj, assay = "RNA", slot = "counts")[, overlap_cells]
  meta <- int_obj@meta.data[overlap_cells, ]
  # Create new Seurat object
  overlap_obj <- CreateSeuratObject(counts, meta.data = meta)
  
  # Perform integration with subset
  subset.list <- c(perserved_obj, overlap_obj)
  subsetint_obj <- seuratv4_integrate(subset.list)
  subsetint_obj <- seurat_clustering(subsetint_obj)
  
  return(subsetint_obj)
}


# Load marker rankings
get_rankings <- function(filepath, top_genes) {
  df <- read.csv(filepath, row.names = 1)
  # Rank by minimum p-value (smaller p → higher rank)
  df$rank <- rank(df$p_val, ties.method = "average")
  # Initialize named vector of NA for all top_genes
  ranks <- setNames(rep(NA, length(top_genes)), top_genes)
  # Fill in only genes that are actually present in the data
  present_genes <- intersect(top_genes, rownames(df))
  ranks[present_genes] <- df[present_genes, "rank"]
  
  return(ranks)
}


# Compare marker gene rankings
compare_rankings <- function(celltype, top_genes,
                             v4_3p.rankings, v4_5p.rankings,
                             subset_3p.rankings, subset_5p.rankings) {
  # Construct dataframe
  rank_comparison <- data.frame(
    gene = top_genes,
    balanced_rank = 1:10,
    v4_3p = v4_3p.rankings[top_genes],
    v4_5p = v4_5p.rankings[top_genes],
    subset_3p = subset_3p.rankings[top_genes],
    subset_5p = subset_5p.rankings[top_genes]
  )
  # Compute Spearman correlations
  spearman_corr <- c(
    cor(rank_comparison$balanced_rank, rank_comparison$v4_3p, method = "spearman"),
    cor(rank_comparison$balanced_rank, rank_comparison$v4_5p, method = "spearman"),
    cor(rank_comparison$balanced_rank, rank_comparison$subset_3p, method = "spearman"),
    cor(rank_comparison$balanced_rank, rank_comparison$subset_5p, method = "spearman")
  )
  # Compute average absolute rank differences
  avg_rank_diff <- c(
    mean(abs(rank_comparison$balanced_rank - rank_comparison$v4_3p)),
    mean(abs(rank_comparison$balanced_rank - rank_comparison$v4_5p)),
    mean(abs(rank_comparison$balanced_rank - rank_comparison$subset_3p)),
    mean(abs(rank_comparison$balanced_rank - rank_comparison$subset_5p))
  )
  
  # Append to dataframe
  metrics_df <- data.frame(
    metric = c("spearman", "avg_rank_diff"),
    v4_3p = c(spearman_corr[1], avg_rank_diff[1]),
    v4_5p = c(spearman_corr[2], avg_rank_diff[2]),
    subset_3p = c(spearman_corr[3], avg_rank_diff[3]),
    subset_5p = c(spearman_corr[4], avg_rank_diff[4])
  )
  
  # Save to CSV
  save_path <- "results/marker_gene_stability"
  compare_name <- paste0(celltype, "_compare")
  metric_name <- paste0(celltype, "_metric")
  write.csv(rank_comparison, 
            file = file.path(save_path, paste0(compare_name, ".csv")), row.names = FALSE)
  write.csv(metrics_df, 
            file = file.path(save_path, paste0(metric_name, ".csv")), row.names = FALSE)
  
  return(list(s = rank_comparison, metrics = metrics_df))
}


# Extract labels and clusters for ARI
ari_prep <- function(annotations, seurat_clusters, celltype, filename) {
  df <- data.frame(true = annotations, predicted = seurat_clusters)
  path <- "results/balanced_metrics/data/"
  write.csv(df, file = paste0(path, filename, "_", celltype, ".csv"), row.names = FALSE)
}


# Get confusion matrix
confusion_matrix <- function(seurat_obj, pred_meta, true_meta = "anno", 
                             celltype, dataset) {
  pred <- seurat_obj[[pred_meta]][, 1]
  true <- seurat_obj[[true_meta]][, 1]
  
  tp <- sum(pred == celltype & true == celltype)
  tn <- sum(pred != celltype & true != celltype)
  fp <- sum(pred == celltype & true != celltype)
  fn <- sum(pred != celltype & true == celltype)
  
  accuracy <- (tp + tn) / (tp + tn + fp + fn)
  precision <- if ((tp + fp) == 0) NA else tp / (tp + fp)
  recall <- if ((tp + fn) == 0) NA else tp / (tp + fn)
  f1 <- if ((2 * tp + fp + fn) == 0) NA else (2 * tp) / (2 * tp + fp + fn)
  
  cat(paste0("==== ", dataset, ": ", celltype,  " ====\n"))
  cat("TP:", tp, "\n")
  cat("TN:", tn, "\n")
  cat("FP:", fp, "\n")
  cat("FN:", fn, "\n")
  cat("Accuracy:", round(accuracy, 4), "\n")
  cat("Precision:", round(precision, 4), "\n")
  cat("Recall:", round(recall, 4), "\n")
  cat("F1:", round(f1, 4), "\n")
}
