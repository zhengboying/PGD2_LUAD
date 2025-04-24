
RunKEGG <- function(FileName = NULL, FilePath = NULL, ResultPath = NULL,
                    UpGene = TRUE, DownGene = TRUE, organism = "hsa",PDF = OUTPUTPDF) {
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
      UpGeneAnnoKEGG <- enrichKEGG(gene = GeneList, keyType = "kegg", organism = organism,
                                   pvalueCutoff = 1, qvalueCutoff = 1,)
      
      if (!is.null(UpGeneAnnoKEGG)) {
        AnnoResult <- as.data.frame(UpGeneAnnoKEGG)
        if (organism == "mmu") {
          AnnoResult$Description <- gsub(" \\- Mus musculus \\(house mouse\\)", "", AnnoResult$Description)
        }
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_UpGene_KEGG.csv"),
                  quote = TRUE, row.names = FALSE)
        
        AnnoResult <- filter(AnnoResult, pvalue < 0.05)
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_UpGene_KEGG_pvalue_0.05.csv"),
                  quote = TRUE, row.names = FALSE)
        
        ## Filter top20 result with the smallest pvalue
        AnnoResult <- as.data.frame(UpGeneAnnoKEGG)
        if (organism == "mmu") {
          AnnoResult$Description <- gsub(" \\- Mus musculus \\(house mouse\\)", "", AnnoResult$Description)
        }
        AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
        
        p1 <- ggplot(AnnoResult, aes(x = reorder(Description, pvalue), y = -log10(pvalue))) &
          geom_bar(fill = "#F8766D", color = "#F8766D", position = "stack", stat = "identity", width = 0.8) &
          labs(x = NULL) & scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none", panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(p1, file = paste0(ResultPath, i, "_UpGene_KEGG_Top20_Barplot.png"),
               width = 6, height = 7)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_UpGene_KEGG_Top20_Barplot.pdf"),
                 width = 6, height = 7)
        }
        
        
        p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -pvalue))) &
          geom_point(aes(size = Count), color = "#F8766D") &
          labs(y = NULL) & scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none",
                panel.grid = element_line(colour = "grey92"), 
                panel.grid.minor = element_line(linewidth = rel(0.5)),
                panel.background = element_rect(fill = "white"))
        ggsave(p1, file = paste0(ResultPath, i, "_UpGene_KEGG_Top20_Dotplot.png"),
               width = 8, height = 6)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_UpGene_KEGG_Top20_Dotplot.pdf"),
                 width = 8, height = 6)
        }
        
        ##
        ##save(UpGeneAnnoKEGG, file = paste0("./Rds/", ResultPath, i, "_UpGene_KEGG.Rda"))
      }
    }
    
    if (DownGene) {
      ## Get down gene
      GeneList <- filter(CalculateResult, avg_log2FC < 0)$ENTREZID
      GeneList <- GeneList[!is.na(GeneList)]
      
      ## Annotation
      DownGeneAnnoKEGG <- enrichKEGG(gene = GeneList, keyType = "kegg", organism = organism,
                                     pvalueCutoff = 1, qvalueCutoff = 1,)
      
      if (!is.null(DownGeneAnnoKEGG)) {
        AnnoResult <- as.data.frame(DownGeneAnnoKEGG)
        if (organism == "mmu") {
          AnnoResult$Description <- gsub(" \\- Mus musculus \\(house mouse\\)", "", AnnoResult$Description)
        }
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_DownGene_KEGG.csv"),
                  quote = TRUE, row.names = FALSE)
        
        AnnoResult <- filter(AnnoResult, pvalue < 0.05)
        AnnoResult <- arrange(AnnoResult, pvalue)
        write.csv(AnnoResult, file = paste0(ResultPath, i, "_DownGene_KEGG_pvalue_0.05.csv"),
                  quote = TRUE, row.names = FALSE)
        
        ## Filter top20 result with the smallest pvalue
        AnnoResult <- as.data.frame(DownGeneAnnoKEGG)
        if (organism == "mmu") {
          AnnoResult$Description <- gsub(" \\- Mus musculus \\(house mouse\\)", "", AnnoResult$Description)
        }
        AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
        
        p1 <- ggplot(AnnoResult, aes(x = reorder(Description, pvalue), y = -log10(pvalue))) &
          geom_bar(fill = "#F8766D", color = "#F8766D", position = "stack", stat = "identity", width = 0.8) &
          labs(x = NULL) & scale_x_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none", panel.background = element_rect(fill = "white"),
                axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
        ggsave(p1, file = paste0(ResultPath, i, "_DownGene_KEGG_Top20_Barplot.png"),
               width = 6, height = 7)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_DownGene_KEGG_Top20_Barplot.pdf"),
                 width = 6, height = 7)
        }
        
        p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -pvalue))) &
          geom_point(aes(size = Count), color = "#F8766D") &
          labs(y = NULL) & scale_y_discrete(labels = function(x) str_wrap(x, width = 50)) &
          PlotTheme$Box & PlotTheme$FontColor &
          theme(legend.position = "none",
                panel.grid = element_line(colour = "grey92"), 
                panel.grid.minor = element_line(linewidth = rel(0.5)),
                panel.background = element_rect(fill = "white"))
        ggsave(p1, file = paste0(ResultPath, i, "_DownGene_KEGG_Top20_Dotplot.png"),
               width = 8, height = 6)
        if (PDF == TRUE){
           ggsave(p1, file = paste0(ResultPath, i, "_DownGene_KEGG_Top20_Dotplot.pdf"),
                  width = 8, height = 6)
        }
        
        ##
        ##save(DownGeneAnnoKEGG, file = paste0("./Rds/", ResultPath, i, "_DownGene_KEGG.Rda"))
      }
    }
  }
}
