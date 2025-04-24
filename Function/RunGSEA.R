
RunGSEA <- function(FileName = NULL, FilePath = NULL, ResultPath = NULL, PDF = FALSE,
                    organism = "Human",
                    GESA_H = TRUE, GESA_C2 = TRUE, GESA_C3 = TRUE,
                    GESA_C5 = TRUE, GESA_C8 = TRUE, 
                    GESA_Custom = FALSE, CustomName = "M1vsM2",CustomFile = "../../Database/gmt/h.all.v2024.1.Hs.symbols.gmt") {
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
    
    ## Sort by log2FC
    GeneList <- as.numeric(CalculateResult$avg_log2FC)
    names(GeneList) <- as.character(CalculateResult$GeneName)
    GeneList <- sort(GeneList, decreasing = TRUE)
    
    ## H
    if (GESA_H) {
      if (organism == "Human") {
        MsigDB_H <- read.gmt("../../Database/gmt/h.all.v2024.1.Hs.symbols.gmt")
      }
      
      if (organism == "Mouse") {
        MsigDB_H <- read.gmt("../../Database/gmt/mh.all.v2024.1.Mm.symbols.gmt")
      }
      
      AnnoH <- GSEA(GeneList, TERM2GENE = MsigDB_H, pvalueCutoff = 1, eps = 0,
                    seed = 16)
      AnnoResult <- as.data.frame(AnnoH)
      write.csv(AnnoResult, file = paste0(ResultPath, i, "_GSEA_H.csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Filter top20 H result with the smallest pvalue
      AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
      AnnoResult <- arrange(AnnoResult, pvalue)
      AnnoResult$Order <- seq(1:nrow(AnnoResult))
      AnnoResult$Count <- (str_count(AnnoResult$core_enrichment, "/") + 1)
      
      p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -Order), color = NES, size = Count)) &
        geom_point() & labs(y = NULL) & scale_colour_gradient(low = "blue", high = "red") &
        PlotTheme$Box & PlotTheme$FontColor &
        theme(panel.grid = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(linewidth = rel(0.5)),
              panel.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white", colour = NA))
      ggsave(p1, file = paste0(ResultPath, i, "_GSEA_H.png"),
             width = 7, height = 6)
      if (PDF == TRUE){
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_H.pdf"),
               width = 7, height = 6)
      }
      
      ## save Rda
      save(AnnoH, file = paste0("./Rds/", ResultPath, i, "_GSEA_H.Rda"))
    }
    
    ## C2
    if (GESA_C2) {
      if (organism == "Human") {
        MsigDB_C2 <- read.gmt("../../Database/gmt/c2.cp.v2024.1.Hs.symbols.gmt")
      }
      
      if (organism == "Mouse") {
        MsigDB_C2 <- read.gmt("../../Database/gmt/m2.cp.v2024.1.Mm.symbols.gmt")
      }
      
      AnnoC2 <- GSEA(GeneList, TERM2GENE = MsigDB_C2, pvalueCutoff = 1, eps = 0,
                     seed = 16)
      AnnoResult <- as.data.frame(AnnoC2)
      write.csv(AnnoResult, file = paste0(ResultPath, i, "_GSEA_C2.csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Filter top20 C2 result with the smallest pvalue
      AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
      AnnoResult <- arrange(AnnoResult, pvalue)
      AnnoResult$Order <- seq(1:nrow(AnnoResult))
      AnnoResult$Count <- (str_count(AnnoResult$core_enrichment, "/") + 1)
      
      p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -Order), color = NES, size = Count)) &
        geom_point() & labs(y = NULL) & scale_colour_gradient(low = "blue", high = "red") &
        PlotTheme$Box & PlotTheme$FontColor &
        theme(panel.grid = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(linewidth = rel(0.5)),
              panel.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white", colour = NA))
      if (max(str_length(AnnoResult$ID)) > 100){
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C2.png"),
               width = max(str_length(AnnoResult$ID))*0.11, height = 6)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C2.pdf"),
                 width = max(str_length(AnnoResult$ID))*0.11, height = 6)
        }
      } else{
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C2.png"),
               width = 9, height = 6)
        if (PDF == TRUE){
            ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C2.pdf"),
                   width = 9, height = 6)
        }
      }
      
      ## save Rda
      save(AnnoC2, file = paste0("./Rds/", ResultPath, i, "_GSEA_C2.Rda"))
    }
    
    ## C3
    if (GESA_C3) {
      if (organism == "Human") {
        MsigDB_C3 <- read.gmt("../../Database/gmt/c3.all.v2024.1.Hs.symbols.gmt")
      }
      
      if (organism == "Mouse") {
        MsigDB_C3 <- read.gmt("../../Database/gmt/m3.all.v2024.1.Mm.symbols.gmt")
      }
      
      AnnoC3 <- GSEA(GeneList, TERM2GENE = MsigDB_C3, pvalueCutoff = 1, eps = 0,
                     seed = 16)
      AnnoResult <- as.data.frame(AnnoC3)
      write.csv(AnnoResult, file = paste0(ResultPath, i, "_GSEA_C3.csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Filter top20 C3 result with the smallest pvalue
      AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
      AnnoResult <- arrange(AnnoResult, pvalue)
      AnnoResult$Order <- seq(1:nrow(AnnoResult))
      AnnoResult$Count <- (str_count(AnnoResult$core_enrichment, "/") + 1)
      
      p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -Order), color = NES, size = Count)) &
        geom_point() & labs(y = NULL) & scale_colour_gradient(low = "blue", high = "red") &
        PlotTheme$Box & PlotTheme$FontColor &
        theme(panel.grid = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(linewidth = rel(0.5)),
              panel.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white", colour = NA))
      if (max(str_length(AnnoResult$ID)) > 100){
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C3.png"),
               width = max(str_length(AnnoResult$ID))*0.11, height = 6)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C3.pdf"),
                 width = max(str_length(AnnoResult$ID))*0.11, height = 6)
        }
      } else{
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C3.png"),
               width = 9, height = 6)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C3.pdf"),
                 width = 9, height = 6)
        }
      }
      
      ## save Rda
      save(AnnoC3, file = paste0("./Rds/", ResultPath, i, "_GSEA_C3.Rda"))
    }
    
    ## C5
    if (GESA_C5) {
      if (organism == "Human") {
        MsigDB_C5 <- read.gmt("../../Database/gmt/c5.all.v2024.1.Hs.symbols.gmt")
      }
      
      if (organism == "Mouse") {
        MsigDB_C5 <- read.gmt("../../Database/gmt/m5.all.v2024.1.Mm.symbols.gmt")
      }
      
      AnnoC5 <- GSEA(GeneList, TERM2GENE = MsigDB_C5, pvalueCutoff = 1, eps = 0,
                     seed = 16)
      AnnoResult <- as.data.frame(AnnoC5)
      write.csv(AnnoResult, file = paste0(ResultPath, i, "_GSEA_C5.csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Filter top20 C5 result with the smallest pvalue
      AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
      AnnoResult <- arrange(AnnoResult, pvalue)
      AnnoResult$Order <- seq(1:nrow(AnnoResult))
      AnnoResult$Count <- (str_count(AnnoResult$core_enrichment, "/") + 1)
      
      p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -Order), color = NES, size = Count)) &
        geom_point() & labs(y = NULL) & scale_colour_gradient(low = "blue", high = "red") &
        PlotTheme$Box & PlotTheme$FontColor &
        theme(panel.grid = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(linewidth = rel(0.5)),
              panel.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white", colour = NA))
      if (max(str_length(AnnoResult$ID)) > 100){
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C5.png"),
               width = max(str_length(AnnoResult$ID))*0.11, height = 6)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C5.pdf"),
                 width = max(str_length(AnnoResult$ID))*0.11, height = 6)
        }
      } else{
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C5.png"),
               width = 9, height = 6)
        if (PDF == TRUE){
          ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C5.pdf"),
                 width = max(str_length(AnnoResult$ID))*0.11, height = 6)
        }
      }
      
      ## save Rda
      save(AnnoC5, file = paste0("./Rds/", ResultPath, i, "_GSEA_C5.Rda"))
    }
    
    ## C8
    if (GESA_C8) {
      if (organism == "Human") {
        MsigDB_C8 <- read.gmt("../../Database/gmt/c8.all.v2024.1.Hs.symbols.gmt")
      }
      
      if (organism == "Mouse") {
        MsigDB_C8 <- read.gmt("../../Database/gmt/m8.all.v2024.1.Mm.symbols.gmt")
      }
      
      AnnoC8 <- GSEA(GeneList, TERM2GENE = MsigDB_C8, pvalueCutoff = 1, eps = 0,
                     seed = 16)
      AnnoResult <- as.data.frame(AnnoC8)
      write.csv(AnnoResult, file = paste0(ResultPath, i, "_GSEA_C8.csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Filter top20 C8 result with the smallest pvalue
      AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
      AnnoResult <- arrange(AnnoResult, pvalue)
      AnnoResult$Order <- seq(1:nrow(AnnoResult))
      AnnoResult$Count <- (str_count(AnnoResult$core_enrichment, "/") + 1)
      
      p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -Order), color = NES, size = Count)) &
        geom_point() & labs(y = NULL) & scale_colour_gradient(low = "blue", high = "red") &
        PlotTheme$Box & PlotTheme$FontColor &
        theme(panel.grid = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(linewidth = rel(0.5)),
              panel.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white", colour = NA))
      ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C8.png"),
             width = 8, height = 6)
      if (PDF == TRUE){
         ggsave(p1, file = paste0(ResultPath, i, "_GSEA_C8.pdf"),
                width = 8, height = 6)
      }
      
      ## save Rda
      save(AnnoC8, file = paste0("./Rds/", ResultPath, i, "_GSEA_C8.Rda"))
    }
    ## Custom gmt
    if (GESA_Custom) {
        MsigDB_Custom <- read.gmt(CustomFile)
      
      
      AnnoCustom <- GSEA(GeneList, TERM2GENE = MsigDB_Custom, pvalueCutoff = 1, eps = 0,
                     seed = 16)
      AnnoResult <- as.data.frame(AnnoCustom)
      write.csv(AnnoResult, file = paste0(ResultPath, i, "_GSEA_", CustomName, ".csv"),
                quote = TRUE, row.names = FALSE)
      
      ## Filter top20 Custom result with the smallest pvalue
      AnnoResult <- top_n(AnnoResult, n = 20, wt = -pvalue)
      AnnoResult <- arrange(AnnoResult, pvalue)
      AnnoResult$Order <- seq(1:nrow(AnnoResult))
      AnnoResult$Count <- (str_count(AnnoResult$core_enrichment, "/") + 1)
      
      p1 <- ggplot(AnnoResult, aes(x = -log10(pvalue), y = reorder(Description, -Order), color = NES, size = Count)) &
        geom_point() & labs(y = NULL) & scale_colour_gradient(low = "blue", high = "red") &
        PlotTheme$Box & PlotTheme$FontColor &
        theme(panel.grid = element_line(colour = "grey92"), 
              panel.grid.minor = element_line(linewidth = rel(0.5)),
              panel.background = element_rect(fill = "white"),
              legend.key = element_rect(fill = "white", colour = NA))
      ggsave(p1, file = paste0(ResultPath, i, "_GSEA_", CustomName, ".png"),
             width = 8, height = 6)
      if (PDF == TRUE){
        ggsave(p1, file = paste0(ResultPath, i, "_GSEA_", CustomName, ".pdf"),
               width = 8, height = 6)
      }
      
      ## save Rda
      save(AnnoCustom, file = paste0("./Rds/", ResultPath, i, "_GSEA_", CustomName, ".Rda"))
    }
    }
    
    
    
  }
}
