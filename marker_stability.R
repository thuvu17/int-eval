# Analyze marker gene stability results

# Load top marker gene names
# Balanced
commonpath <- "results/marker_gene_stability"
filepath <- paste0(commonpath, "/balanced/", celltype.keep, "_balanced", ".csv")
balanced_top10 <- rownames(read.csv(filepath, row.names = 1))[1:10]
# Normal Seurat v4
filepath <- paste0(commonpath, "/normal_v4/", celltype.keep, "_normal_v4_3p.csv")
v4_top10.3p <- rownames(read.csv(filepath, row.names = 1))[1:10]
filepath <- paste0(commonpath, "/normal_v4/", celltype.keep, "_normal_v4_5p.csv")
v4_top10.5p <- rownames(read.csv(filepath, row.names = 1))[1:10]
# Subset
filepath <- paste0(commonpath, "/subset/", celltype.keep, "_subset_3p.csv")
subset_top10.3p <- rownames(read.csv(filepath, row.names = 1))[1:10]
filepath <- paste0(commonpath, "/subset/", celltype.keep, "_subset_5p.csv")
subset_top10.5p <- rownames(read.csv(filepath, row.names = 1))[1:10]

df_compare <- data.frame(
  gene_rank = 1:10,
  balanced = balanced_top10,
  normal_v4_3p = v4_top10.3p,
  normal_v4_5p = v4_top10.5p,
  subset_3p = subset_top10.3p,
  subset_5p = subset_top10.5p
)
