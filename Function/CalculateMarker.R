# Define a function for marker analysis
CalculateMarker <- function(object = scObject, assay.use = DefaultAssay(scObject), slot.use = "data",
                            ClusterName = NULL, ClusterName.1 = NULL,
                            min.pct = 0.01, logfc.threshold = 0.01,
                            max.cells.per.ident = Inf, only.pos = TRUE,
                            test.use = "wilcox",
                            min.cells.feature = 3, min.cells.group = 3,
                            annotation = AnnoList) {
  # Check if object is a Seurat object
  stopifnot(is(object, "Seurat"))
  
  # Check if ClusterName and ClusterName.1 are provided
  stopifnot(!is.null(ClusterName), !is.null(ClusterName.1))
  
  # Set the identity of cells in the Seurat object to ClusterName
  Idents(object) <- ClusterName
  
  # Find marker
  AllMarker <- FindMarkers(object, assay = assay.use, slot = slot.use,
                           ident.1 = ClusterName.1, min.pct = min.pct, logfc.threshold = logfc.threshold,
                           max.cells.per.ident = max.cells.per.ident, only.pos = only.pos,
                           test.use = test.use, random.seed = 16,
                           min.cells.feature = min.cells.feature, min.cells.group = min.cells.group)
  AllMarker <- cbind(GeneName = row.names(AllMarker), AllMarker)
  
  # Calculate mean expression of CellName.1 and CellName.2（GetAssayData>LayerData）
  Base_Mean1 <- as.data.frame(unlist(rowMeans(expm1(LayerData(subset(object, idents = ClusterName.1, return.null = TRUE), assay = "RNA", layer = "data")))))
  Base_Mean1 <- cbind(GeneName = row.names(Base_Mean1), Base_Mean1)
  Base_Mean2 <- as.data.frame(unlist(rowMeans(expm1(LayerData(subset(object, idents = ClusterName.1, return.null = TRUE, invert = TRUE), assay = "RNA", layer = "data")))))
  Base_Mean2 <- cbind(GeneName = row.names(Base_Mean2), Base_Mean2)
  
  # Combine Base_Mean1 and Base_Mean2 into a single data frame
  Base_Mean1 <- merge(x = Base_Mean1, y = Base_Mean2, by = "GeneName", all = TRUE)
  colnames(Base_Mean1) <- c("GeneName", "exp.1", "exp.2")
  
  # Merge AllMarker with Base_Mean
  AllMarker <- merge(x = AllMarker, y = Base_Mean1, by = "GeneName", all.x = TRUE)
  
  if (!is.null(annotation)) {
    # Merge AllMarker with annotation
    AllMarker <- merge(x = AllMarker, y = annotation, by.x = "GeneName", by.y = "Symbol", all.x = TRUE)
    
    # Reorder columns in AllMarker
    AllMarker <- AllMarker[, c(1, 9, 10, 7, 8, 3, 2, 6, 4, 5, 11)]
  } 
  
  # Sort AllMarker by absolute average log2 fold change
  AllMarker <- AllMarker[order(abs(AllMarker$avg_log2FC), decreasing = TRUE), ]
}
