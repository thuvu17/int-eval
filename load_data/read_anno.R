# Cell type annotation was provided for 3' PBMC
# Read in metadata: batch, cell type = annotation, barcode from
# https://doi.org/10.1093/bioinformatics/btz625 (Polanski et.al.)
meta <- read.csv("data/polanski_metadata.csv", stringsAsFactors = FALSE)
head(meta)

# Subset metadata to batch 3'
meta_3p <- subset(meta, batch == "0")

# Add cell type to matching barcodes
meta_3p$barcode <- sub("-0$", "-1", meta_3p$barcode)
pbmc_3p$celltype <- meta_3p$Cell.type[match(colnames(pbmc_3p), meta_3p$barcode)]

# 8090 cells match, 283 mismatch. Exclude the ones without matching cell type
pbmc_3p <- subset(pbmc_3p, subset = !is.na(celltype))
remove(meta_3p)

# ==============================================================================
# Annotate 5' PBMC with Azimuth
# https://azimuth.hubmapconsortium.org/
meta_5p <- read.delim("data/azimuth_pred.tsv", stringsAsFactors = FALSE)
meta_5p$cell <- paste0("data_5p-", meta_5p$cell)

# Only keep the matching barcodes
common_barcodes <- intersect(meta_5p$cell, colnames(pbmc_5p))
meta_5p <- meta_5p[meta_5p$cell %in% common_barcodes, ]
pbmc_5p <- subset(pbmc_5p, cells = common_barcodes)

# Translating Azimuth annotation to match 3'
l1anno_map <- c(
  "B" = "B cell",
  "CD4 T" = "CD4 T cell",
  "CD8 T" = "CD8 T cell",
  "NK" = "NK cell"
)

l2anno_map <- c(
  "pDC" = "Plasmacytoid dendritic cell",
  "HSPC" = "Hematopoietic stem cell",
  "CD14 Mono" = "Monocyte_CD14",
  "CD16 Mono" = "Monocyte_FCGR3A",
  "Platelet" = "Megakaryocyte"
)

meta_5p$celltype <- l1anno_map[meta_5p$predicted.celltype.l1]
# Use l2 mapping to override NA values
meta_5p$celltype[is.na(meta_5p$celltype)] <- l2anno_map[meta_5p$predicted.celltype.l2[is.na(meta_5p$celltype)]]

# 999 doesn't have the matching cell types so exclude them
meta_5p <- meta_5p[!is.na(meta_5p$celltype), ]

# Only keep barcodes in both meta_5p and pbmc_5p
common_barcodes <- intersect(meta_5p$cell, colnames(pbmc_5p))
meta_5p <- meta_5p[meta_5p$cell %in% common_barcodes, ]
pbmc_5p <- subset(pbmc_5p, cells = common_barcodes)
meta_5p <- meta_5p[match(colnames(pbmc_5p), meta_5p$cell), ]

# Translate over cell type annotation
pbmc_5p$celltype <- meta_5p$celltype
remove(l1anno_map, l2anno_map, common_barcodes)
