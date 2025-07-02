# Read in metadata: batch, cell type = annotation, barcode
meta <- read.csv("data/pbmc_metadata.csv", stringsAsFactors = FALSE)
head(meta)

# Subset metadata to batch 3'
meta_3p <- subset(meta, batch == "0")

# Add cell type to matching barcodes
meta_3p$barcode <- sub("-0$", "-1", meta_3p$barcode)
pbmc_3p$celltype <- meta_3p$Cell.type[match(colnames(pbmc_3p), meta_3p$barcode)]

# 8090 cells match, 283 mismatch. Exclude the ones without matching cell type
pbmc_3p <- subset(pbmc_3p, subset = !is.na(celltype))
