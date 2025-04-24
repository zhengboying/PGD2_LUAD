
ClusterStat <- function(object = scObject, percent = FALSE, RowType = "Celltype", ColType = "Sample") {
  RowMeta <- object@meta.data[[RowType]]
  ColMeta <- object@meta.data[[ColType]]
  
  if (percent == FALSE) {
    ClusterStatResult <- as.data.frame.matrix(table(RowMeta, ColMeta))
    
    Total <- as.data.frame(t(colSums(ClusterStatResult)))
    rownames(Total) <- "Total"
    colnames(Total) <- colnames(ClusterStatResult)
    
    ClusterStatResult <- rbind(ClusterStatResult, Total)
    ClusterStatResult <- cbind(rownames(ClusterStatResult), ClusterStatResult)
    colnames(ClusterStatResult)[1] <- RowType
  } else {
    ClusterStatResult <- as.data.frame.matrix(prop.table(table(RowMeta, ColMeta), 2))
    
    Total <- as.data.frame(t(rep(1, ncol(ClusterStatResult))))
    rownames(Total) <- "Total"
    colnames(Total) <- colnames(ClusterStatResult)
    
    ClusterStatResult <- rbind(ClusterStatResult, Total)
    for (i in colnames(ClusterStatResult)) {
      ClusterStatResult[, i] <- scales::label_percent(accuracy = 0.01)(ClusterStatResult[, i])
    }
    
    ClusterStatResult <- cbind(rownames(ClusterStatResult), ClusterStatResult)
    colnames(ClusterStatResult)[1] <- RowType
  }
  
  return(ClusterStatResult)
}