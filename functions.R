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
