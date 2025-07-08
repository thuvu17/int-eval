# Analyze marker gene stability results


# Load top marker gene names
commonpath <- "results/marker_gene_stability"
filepath <- paste0(commonpath, "/balanced/", celltype.keep, "_balanced", ".csv")
balanced_df <- read.csv(filepath, row.names = 1)
balanced_df$rank <- rank(balanced_df$minimump_p_val, ties.method = "average")
balanced_df_sorted <- balanced_df[order(balanced_df$rank), ]
balanced_df_top10 <- rownames(balanced_df_sorted)[1:10]

# Get top marker gene rankings in v4 and subset
# Normal Seurat v4
filepath <- paste0(commonpath, "/normal_v4/", celltype.keep, "_normal_v4_3p.csv")
v4_3p.rankings <- get_rankings(filepath, balanced_df_top10)
filepath <- paste0(commonpath, "/normal_v4/", celltype.keep, "_normal_v4_5p.csv")
v4_5p.rankings <- get_rankings(filepath, balanced_df_top10)
# Subset
filepath <- paste0(commonpath, "/subset/", celltype.keep, "_subset_3p.csv")
subset_3p.rankings <- get_rankings(filepath, balanced_df_top10)
filepath <- paste0(commonpath, "/subset/", celltype.keep, "_subset_5p.csv")
subset_5p.rankings <- get_rankings(filepath, balanced_df_top10)
remove(filepath)

# Compare rankings between balanced and perturbed
rank_comparison <- data.frame(
  gene = balanced_df_top10,
  balanced_rank = 1:10,
  v4_3p = v4_3p.rankings[balanced_df_top10],
  v4_5p = v4_5p.rankings[balanced_df_top10],
  subset_3p = subset_3p.rankings[balanced_df_top10],
  subset_5p = subset_5p.rankings[balanced_df_top10]
)
