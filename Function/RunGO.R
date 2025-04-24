
RunGO <- function(FileName = NULL, FilePath = NULL, ResultPath = NULL,PDF = FALSE,
                  UpGene = TRUE, DownGene = TRUE, OrgDb = org.Hs.eg.db) {
  if (is.null(FileName)) {
    stop('argument "FileName" is missing.')
  }
  
  if (is.null(FilePath)) {
    stop('argument "FilePath" is missing.')
  }
  
  if (is.null(ResultPath)) {
    stop('argument "ResultPath" is missing.')
  }
  
  for (i in FileName) {
    ## Load calculate result
    CalculateResult <- read.csv(file = paste0(FilePath, "/", i, ".csv"),
                                header = TRUE)
    
    if (UpGene) {
      ## Get up gene
      GeneList <- filter(CalculateResult, avg_log2FC > 0)$ENTREZID
      GeneList <- GeneList[!is.na(GeneList)]
      
      ## Annotation
      UpGeneAnnoGO <- enrichGO(gene = GeneList, keyType = "ENTREZID", OrgDb = OrgDb,
                               ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1,
                               readable = TRUE)
      
      if (!is.null(UpGeneAnnoGO)) {
        AnnoResult <- as.data.frame(UpGeneAnnoGO)
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_UpGene_GO.csv"),
                  quote = TRUE, row.names = FALSE)
        
        AnnoResult <- filter(AnnoResult, pvalue < 0.05)
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_UpGene_GO_pvalue_0.05.csv"),
                  quote = TRUE, row.names = FALSE)
        
        ## Filter top10 result with the smallest pvalue
        AnnoResult <- as.data.frame(UpGeneAnnoGO)
        AnnoResult <- top_n(group_by(AnnoResult, ONTOLOGY), n = 10, wt = -pvalue)
        
        AnnoResult$ONTOLOGY <- factor(AnnoResult$ONTOLOGY, levels = c("BP", "CC", "MF"))
        p1 <- ggplot(AnnoResult, aes(x = reorder(Description, pvalue), y = -log10(pvalue))) &
          geom_bar(aes(fill = ONTOLOGY, color = ONTOLOGY), position = "stack", stat = "identity", width = 0.8) &
          facet_grid(. ~ ONTOLOGY, scales = "free_x", as.table = FALSE, space = "free") &
          labs(x = NULL) & scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none", panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(p1, file = paste0(ResultPath, i, "_UpGene_GO_Top10_Barplot.png"),
               width = 9, height = 7)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_UpGene_GO_Top10_Barplot.pdf"),
                 width = 9, height = 7)
        }
        
        AnnoResult$ONTOLOGY <- factor(AnnoResult$ONTOLOGY, levels = c("MF", "CC", "BP"))
        p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -pvalue))) &
          geom_point(aes(size = Count, fill = ONTOLOGY, color = ONTOLOGY)) &
          facet_grid(ONTOLOGY ~ ., scales = "free_y", as.table = F, space = "free") &
          labs(y = NULL) & scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none",
                panel.grid = element_line(colour = "grey92"), 
                panel.grid.minor = element_line(linewidth = rel(0.5)),
                panel.background = element_rect(fill = "white"))
        ggsave(p1, file = paste0(ResultPath, i, "_UpGene_GO_Top10_Dotplot.png"),
               width = 8, height = 9)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_UpGene_GO_Top10_Dotplot.pdf"),
                 width = 8, height = 9)
        }
        
        ##
        ##save(UpGeneAnnoGO, file = paste0("./Rds/", ResultPath, i, "_UpGene_GO.Rda"))
      }
    }
    
    if (DownGene) {
      ## Get down gene
      GeneList <- filter(CalculateResult, avg_log2FC < 0)$ENTREZID
      GeneList <- GeneList[!is.na(GeneList)]
      
      ## Annotation
      DownGeneAnnoGO <- enrichGO(gene = GeneList, keyType = "ENTREZID", OrgDb = OrgDb,
                                 ont = "ALL", pvalueCutoff = 1, qvalueCutoff = 1,
                                 readable = TRUE)
      
      if (!is.null(DownGeneAnnoGO)) {
        AnnoResult <- as.data.frame(DownGeneAnnoGO)
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_DownGene_GO.csv"),
                  quote = TRUE, row.names = FALSE)
        
        AnnoResult <- filter(AnnoResult, pvalue < 0.05)
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_DownGene_GO_pvalue_0.05.csv"),
                  quote = TRUE, row.names = FALSE)
        
        ## Filter top10 result with the smallest pvalue
        AnnoResult <- as.data.frame(DownGeneAnnoGO)
        AnnoResult <- top_n(group_by(AnnoResult, ONTOLOGY), n = 10, wt = -pvalue)
        
        AnnoResult$ONTOLOGY <- factor(AnnoResult$ONTOLOGY, levels = c("BP", "CC", "MF"))
        p1 <- ggplot(AnnoResult, aes(x = reorder(Description, pvalue), y = -log10(pvalue))) &
          geom_bar(aes(fill = ONTOLOGY, color = ONTOLOGY), position = "stack", stat = "identity", width = 0.8) &
          facet_grid(. ~ ONTOLOGY, scales = "free_x", as.table = FALSE, space = "free") &
          labs(x = NULL) & scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none", panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(p1, file = paste0(ResultPath, i, "_DownGene_GO_Top10_Barplot.png"),
               width = 9, height = 7)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_DownGene_GO_Top10_Barplot.pdf"),
                 width = 9, height = 7)
        }
        
        AnnoResult$ONTOLOGY <- factor(AnnoResult$ONTOLOGY, levels = c("MF", "CC", "BP"))
        p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -pvalue))) &
          geom_point(aes(size = Count, fill = ONTOLOGY, color = ONTOLOGY)) &
          facet_grid(ONTOLOGY ~ ., scales = "free_y", as.table = F, space = "free") &
          labs(y = NULL) & scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none",
                panel.grid = element_line(colour = "grey92"), 
                panel.grid.minor = element_line(linewidth = rel(0.5)),
                panel.background = element_rect(fill = "white"))
        ggsave(p1, file = paste0(ResultPath, i, "_DownGene_GO_Top10_Dotplot.png"),
               width = 8, height = 9)
        if (PDF == TRUE){
           ggsave(p1, file = paste0(ResultPath, i, "_DownGene_GO_Top10_Dotplot.pdf"),
                  width = 8, height = 9)
        }
        
        ##
        ##save(DownGeneAnnoGO, file = paste0("./Rds/", ResultPath, i, "_DownGene_GO.Rda"))
      }
    }
  }
}
