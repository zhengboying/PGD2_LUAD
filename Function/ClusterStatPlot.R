
ClusterStatPlot <- function(PlotData = ClusterStatResult,
                            RowType = NULL, ColType = NULL, PlotPrefixName = NULL,
                            RowTypeAngle = FALSE, ColTypeAngle = FALSE, AddWidth = 2, 
                            ResultPath = "./",PDF = FALSE) {
  if (is.null(RowType)) {
    stop('argument "RowType" is missing.')
  }
  
  if (is.null(ColType)) {
    stop('argument "ColType" is missing.')
  }
  
  if (is.null(PlotPrefixName)) {
    stop('argument "PlotPrefixName" is missing.')
  }
  
  ##
  PlotData <-ClusterStatResult
  PlotData <- PlotData[-nrow(PlotData),]
  colnames(PlotData)[1] <- "Cluster"
  PlotData <- melt(PlotData, id.vars = "Cluster")
  PlotData$Cluster <- factor(PlotData$Cluster)
  nlevels <- length(levels(PlotData$Cluster))
  #PlotData$Cluster <- factor(PlotData$Cluster, levels = levels(scObject@meta.data[[RowType]]))
  
  ## Plot
  {
    p1 <- ggplot(PlotData, aes(x = Cluster, y = as.numeric(value))) &
      geom_bar(aes(fill = variable), position = "stack", stat = "identity", width = 0.9) &
      labs(x = NULL, y = "Number of Cells") &
      scale_fill_manual(values = PlotTheme$Color) & PlotTheme$Box & PlotTheme$FontColor &
      theme(legend.title = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"))
    if (RowTypeAngle) { p1 <- p1 & theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) }
    ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_Cluster_", ColType, "_Stat.png"),
           bg = "white", width = (0.4*nlevels(p1[["data"]][["Cluster"]]) + AddWidth), height = 6)
    if (PDF == TRUE){
      ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_Cluster_", ColType, "_Stat.pdf"),
             bg = "white", width = (0.4*nlevels(p1[["data"]][["Cluster"]]) + AddWidth), height = 6)
    }
    
    
    p1 <- ggplot(PlotData, aes(x = variable, y = as.numeric(value))) &
      geom_bar(aes(fill = Cluster), position = "stack", stat = "identity", width = 0.9) &
      labs(x = NULL, y = "Number of Cells") &
      scale_fill_manual(values = PlotTheme$Color) & PlotTheme$Box & PlotTheme$FontColor &
      theme(legend.title = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"))
    if (ColTypeAngle) { p1 <- p1 & theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) }
    ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_", ColType, "_Stat.png"),
           bg = "white", width = (0.5*nlevels(p1[["data"]][["variable"]]) + 2.5), height = 6)
    if (PDF == TRUE){
      ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_", ColType, "_Stat.pdf"),
             bg = "white", width = (0.5*nlevels(p1[["data"]][["variable"]]) + 2.5), height = 6)
    }
  }
  
  if (length(unique(PlotData$variable)) > 1) {
    p1 <- ggplot(PlotData, aes(x = Cluster, y = as.numeric(value))) &
      geom_bar(aes(fill = variable), position = "fill", stat = "identity", width = 0.9) &
      labs(x = NULL, y = "Fraction of Cells") &
      scale_fill_manual(values = PlotTheme$Color) & PlotTheme$Box & PlotTheme$FontColor &
      theme(legend.title = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"))
    if (RowTypeAngle) { p1 <- p1 & theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) }
    ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_Cluster_", ColType, "_Stat_Percent.png"),
           bg = "white", width = (0.4*nlevels(p1[["data"]][["Cluster"]]) + AddWidth), height = 6)
    if (PDF == TRUE){
      ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_Cluster_", ColType, "_Stat_Percent.pdf"),
             bg = "white", width = (0.4*nlevels(p1[["data"]][["Cluster"]]) + AddWidth), height = 6)
    }
    
    p1 <- ggplot(PlotData, aes(x = variable, y = as.numeric(value))) &
      geom_bar(aes(fill = Cluster), position = "fill", stat = "identity", width = 0.9) &
      labs(x = NULL, y = "Fraction of Cells") &
      scale_fill_manual(values = PlotTheme$Color) & PlotTheme$Box & PlotTheme$FontColor &
      theme(legend.title = element_blank(), panel.grid = element_blank(),
            panel.background = element_rect(fill = "white"))
    if (ColTypeAngle) { p1 <- p1 & theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) }
    ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_", ColType, "_Stat_Percent.png"),
           bg = "white", width = (0.5*nlevels(p1[["data"]][["variable"]]) + 2.5), height = 6)
    if (PDF == TRUE){
      ggsave(p1, file = paste0(ResultPath, PlotPrefixName, "_", ColType, "_Stat_Percent.pdf"),
             bg = "white", width = (0.5*nlevels(p1[["data"]][["variable"]]) + 2.5), height = 6)
    }
  }
}