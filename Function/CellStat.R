

## CellStat
CellStat <- function(object = scObject, Type = "orig.ident", assay.use = NULL) {
  if (!is(object, "Seurat")) {
    stop("object must be a Seurat object.")
  }
  
  if (is.null(assay.use)) {
    stop('assay.use has to one form the list: "VisiumHD", "Spatial", "RNA", "ATAC", "ARC".')
  }
  
  if (!((assay.use == "VisiumHD") | (assay.use == "Spatial") | (assay.use == "RNA") | (assay.use == "ATAC") | (assay.use == "ARC"))) {
    stop('assay.use has to one form the list: "VisiumHD", "Spatial", "RNA", "ATAC", "ARC".')
  }
  
  Idents(object) <- Type
  
  ## ARC
  if (assay.use == "ARC") {
    ## Create an empty matrix
    df <- data.frame(matrix(nrow = nlevels(object), ncol = 10))
    colnames(df) <- c("Sample", "CellNumber",
                      "MeanUMI_per_cell", "MedianUMI_per_cell",
                      "MeanGene_per_cell", "MedianGene_per_cell",
                      "MeanRead_per_cell", "MedianRead_per_cell",
                      "MeanPeak_per_cell", "MedianPeak_per_cell")
    df$Sample <- levels(object)
    
    ##
    for(i in 1:nlevels(object)) {
      cell_only = (Idents(object) == levels(object)[i])
      
      df[i, 2] = round(sum(cell_only))
      df[i, 3] = round(mean(object$nCount_RNA[cell_only], na.rm = TRUE))
      df[i, 4] = round(median(object$nCount_RNA[cell_only], na.rm = TRUE))
      df[i, 5] = round(mean(object$nFeature_RNA[cell_only], na.rm = TRUE))
      df[i, 6] = round(median(object$nFeature_RNA[cell_only], na.rm = TRUE))
      df[i, 7] = round(mean(object$nCount_ATAC[cell_only], na.rm = TRUE))
      df[i, 8] = round(median(object$nCount_ATAC[cell_only], na.rm = TRUE))
      df[i, 9] = round(mean(object$nFeature_ATAC[cell_only], na.rm = TRUE))
      df[i, 10] = round(median(object$nFeature_ATAC[cell_only], na.rm = TRUE))
    }
    
    if (nlevels(object) > 1) {
      df <- rbind(df, c("Total",
                        round(length(names(Idents(object)))),
                        round(mean(object$nCount_RNA, na.rm = TRUE)),
                        round(median(object$nCount_RNA, na.rm = TRUE)),
                        round(mean(object$nFeature_RNA, na.rm = TRUE)),
                        round(median(object$nFeature_RNA, na.rm = TRUE)),
                        round(mean(object$nCount_ATAC, na.rm = TRUE)),
                        round(median(object$nCount_ATAC, na.rm = TRUE)),
                        round(mean(object$nFeature_ATAC, na.rm = TRUE)),
                        round(median(object$nFeature_ATAC, na.rm = TRUE))))
    }
  } else {
    ## VisiumHD
    if (assay.use == "VisiumHD") {
      ## Create an empty matrix
      df <- data.frame(matrix(nrow = nlevels(object), ncol = 6))
      colnames(df) <- c("Sample", "BinNumber",
                        "MeanUMI_per_bin", "MedianUMI_per_bin",
                        "MeanGene_per_bin", "MedianGene_per_bin")
      df$Sample <- levels(object)
    }
    
    ## Spatial
    if (assay.use == "Spatial") {
      ## Create an empty matrix
      df <- data.frame(matrix(nrow = nlevels(object), ncol = 6))
      colnames(df) <- c("Sample", "SpotNumber",
                        "MeanUMI_per_spot", "MedianUMI_per_spot",
                        "MeanGene_per_spot", "MedianGene_per_spot")
      df$Sample <- levels(object)
    }
    
    ## RNA
    if (assay.use == "RNA") {
      ## Create an empty matrix
      df <- data.frame(matrix(nrow = nlevels(object), ncol = 6))
      colnames(df) <- c("Sample", "CellNumber",
                        "MeanUMI_per_cell", "MedianUMI_per_cell",
                        "MeanGene_per_cell", "MedianGene_per_cell")
      df$Sample <- levels(object)
    }
    
    ## ATAC
    if (assay.use == "ATAC") {
      ## Create an empty matrix
      df <- data.frame(matrix(nrow = nlevels(object), ncol = 6))
      colnames(df) <- c("Sample", "CellNumber",
                        "MeanRead_per_cell", "MedianRead_per_cell",
                        "MeanPeak_per_cell", "MedianPeak_per_cell")
      df$Sample <- levels(object)
    }
    
    ##
    for(i in 1:nlevels(object)) {
      cell_only = (Idents(object) == levels(object)[i])
      
      df[i, 2] = round(sum(cell_only))
      df[i, 3] = round(mean(object@meta.data[[paste0("nCount_", assay.use)]][cell_only], na.rm = TRUE))
      df[i, 4] = round(median(object@meta.data[[paste0("nCount_", assay.use)]][cell_only], na.rm = TRUE))
      df[i, 5] = round(mean(object@meta.data[[paste0("nFeature_", assay.use)]][cell_only], na.rm = TRUE))
      df[i, 6] = round(median(object@meta.data[[paste0("nFeature_", assay.use)]][cell_only], na.rm = TRUE))
    }
    
    if (nlevels(object) > 1) {
      df <- rbind(df, c("Total",
                        round(length(names(Idents(object)))),
                        round(mean(object@meta.data[[paste0("nCount_", assay.use)]], na.rm = TRUE)),
                        round(median(object@meta.data[[paste0("nCount_", assay.use)]], na.rm = TRUE)),
                        round(mean(object@meta.data[[paste0("nFeature_", assay.use)]], na.rm = TRUE)),
                        round(median(object@meta.data[[paste0("nFeature_", assay.use)]], na.rm = TRUE))))
    }
  }
  
  return(df)
}
