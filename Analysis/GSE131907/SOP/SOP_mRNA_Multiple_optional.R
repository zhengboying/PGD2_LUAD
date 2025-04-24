

## 08.2 res.Marker at -----
suppressMessages({
  library(dplyr)
  library(stringr)
})

source("../../Function/ResChoose.R")
source("../../Function/CalculateMarker.R")

KeyWord = "Harmony"#"RAW"
VarNumber = 15
ReductionNumber = c(0.1, 0.2, 0.4, 0.6)


ResultPath <- c("./4.res.Marker/")

if (str_to_title(SampleSpecies) == "Human") {
  AnnoList <- read.table("../../Database/Human_Ensembl_Symbol_ENTREZID.xls", header = TRUE, sep = "\t", quote = "")
} else if (str_to_title(SampleSpecies) == "Mouse") {
  AnnoList <- read.table("../../Database/Mouse_Ensembl_Symbol_ENTREZID.xls", header = TRUE, sep = "\t", quote = "")
} else {
  stop("SampleSpecies must be one of the following: Human, Mouse.")
}

##CalculateMarker
{
  ## Create result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  
  ##
  ResNubmer <- ResChoose(object = scObject, KeyWord = KeyWord, VarNumber = VarNumber,
                         ReductionNumber = ReductionNumber)
  
  
  TopGene <- c()
  TopNumber <- 10
  
  Idents(scObject) <- paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")
  for (i in levels(scObject)) {
    CalculateResult <- CalculateMarker(object = scObject,
                                       ClusterName = paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony"),
                                       ClusterName.1 = i,
                                       min.pct = 0.01, logfc.threshold = 0.01,
                                       max.cells.per.ident = 3000, only.pos = TRUE,
                                       annotation = AnnoList)
    write.csv(CalculateResult, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Marker.csv"),
              quote = TRUE, row.names = FALSE)
    
    TopGene <- c(TopGene, top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName)
    
    {
      p1 <- VlnPlot(scObject, features = top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName,
                    pt.size = 0.01, assay = "RNA", raster = FALSE) &
        NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor
      ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_ViolinPlot.png"),
             bg = "white", width = 17, height = 10)
      if (OUTPUTPDF){
        ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_ViolinPlot.pdf"),
               bg = "white", width = 17, height = 10)
      }
      
      p1 <- FeaturePlot(scObject, features = top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName, reduction = "umap",
                        raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),
                        pt.size = 0.05, max.cutoff = "q99") &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_FeaturePlot.png"),
             bg = "white", width = 16, height = 10)
      if (OUTPUTPDF){
        ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_FeaturePlot.pdf"),
               bg = "white", width = 16, height = 10)
      }
    }
    
    CalculateResult <- filter(CalculateResult, p_val < 0.05)
    CalculateResult <- filter(CalculateResult, avg_log2FC >= 0.26)
    write.csv(CalculateResult, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Marker_pvalue_0.05_log2FC_0.26.csv"),
              quote = TRUE, row.names = FALSE)
  }
}
##PLOT
{
  {
    p1 <- DotPlot(scObject, features = unique(TopGene), cols = c("white", "red"),
                  cluster.idents = FALSE, assay = "RNA", dot.scale = 4, col.min = 0) &
      labs(x = NULL, y = NULL) & PlotTheme$Box & PlotTheme$FontColor &
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    p1[["data"]][["id"]] <- factor(p1[["data"]][["id"]], levels = rev(levels(p1[["data"]][["id"]])))
    ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Marker_Top", TopNumber, "_Dotplot.png"),
           bg = "white", width = 24, height = 6)
    if (OUTPUTPDF){
      ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Marker_Top", TopNumber, "_Dotplot.pdf"),
             bg = "white", width = 24, height = 6)
    }
    
  }
  
  {
    cluster.averages <- as.data.frame(t(AverageExpression(scObject)[["RNA"]][unique(TopGene),]))
    
    p1 <- pheatmap::pheatmap(cluster.averages, scale = "column", cluster_rows = FALSE, cluster_cols = FALSE,
                             border_color = NA, angle_col = 90)
    ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Marker_Top", TopNumber, "_Heatmap.png"),
           bg = "white", width = 24, height = 6)
    if (OUTPUTPDF){
      ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Marker_Top", TopNumber, "_Heatmap.pdf"),
             bg = "white", width = 24, height = 6)
    }
  }
}




## res.Marker_GO_KEGG -----
suppressMessages({
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
})

source("../../Function/RunGO.R")
source("../../Function/RunKEGG.R")

CalculatePath <- c("./4.res.Marker/")
ResultPathGO <- c("./5.Enrichment/1.GO/")
ResultPathKEGG <- c("./5.Enrichment/2.KEGG/")

##
{
  ## Get the file name containing only significant information
  FileName <- list.files(CalculatePath)
  FileName <- FileName[str_detect(FileName, ".csv")]
  FileName <- FileName[str_detect(FileName, "pvalue", negate = FALSE)]
  FileName <- gsub(".csv", "", FileName)
  
  ## Create GO/KEGG result path
  dir.create(ResultPathGO, showWarnings = FALSE, recursive = TRUE)
  dir.create(ResultPathKEGG, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./Rds/", ResultPathGO), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./Rds/", ResultPathKEGG), showWarnings = FALSE, recursive = TRUE)
  
  if (str_to_title(SampleSpecies) == "Human") {
    suppressMessages({
      library(org.Hs.eg.db)
    })
    
    NeedOrgDb <- org.Hs.eg.db
    NeedOrganism <- "hsa"
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    suppressMessages({
      library(org.Mm.eg.db)
    })
    
    NeedOrgDb <- org.Mm.eg.db
    NeedOrganism <- "mmu"
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  ## Run GO
  RunGO(FileName = FileName,
        PDF = OUTPUTPDF,
        FilePath = CalculatePath,
        ResultPath = ResultPathGO,
        OrgDb = NeedOrgDb,
        UpGene = TRUE, DownGene = FALSE)
  
  ## Run KEGG
  RunKEGG(FileName = FileName,
          PDF = OUTPUTPDF,
          FilePath = CalculatePath,
          ResultPath = ResultPathKEGG,
          organism = NeedOrganism,
          UpGene = TRUE, DownGene = FALSE)
}

## res.Marker_GSEA -----
suppressMessages({
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
})

source("../../Function/RunGSEA.R")

CalculatePath <- c("./4.res.Marker/")
ResultPath <- c("./5.Enrichment/3.GSEA/")

##
{
  ## Get the file name containing all information
  FileName <- list.files(CalculatePath)
  FileName <- FileName[str_detect(FileName, ".csv")]
  FileName <- FileName[str_detect(FileName, "pvalue", negate = TRUE)]
  FileName <- gsub(".csv", "", FileName)
  
  ## Create GESA result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./Rds/", ResultPath), showWarnings = FALSE, recursive = TRUE)
  
  if (str_to_title(SampleSpecies) == "Human" | str_to_title(SampleSpecies) == "Mouse") {
    ## Run GESA
    RunGSEA(FileName = FileName,
            PDF = OUTPUTPDF,
            FilePath = CalculatePath,
            ResultPath = ResultPath,
            organism = str_to_title(SampleSpecies),
            GESA_H = TRUE, GESA_C2 = TRUE, GESA_C3 = TRUE,
            GESA_C5 = TRUE, GESA_C8 = TRUE,
            GESA_Custom = FALSE, CustomName = "M1vsM2", CustomFile = "../../Database/gmt/h.all.v2024.1.Hs.symbols.gmt")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
}


if (str_to_title(SampleSpecies) == "Human") {
  AnnoList <- read.table("/home/kaicheng/file/Human_Ensembl_Symbol_ENTREZID.xls", header = TRUE, sep = "\t", quote = "")
} else if (str_to_title(SampleSpecies) == "Mouse") {
  AnnoList <- read.table("/home/kaicheng/file/Mouse_Ensembl_Symbol_ENTREZID.xls", header = TRUE, sep = "\t", quote = "")
} else {
  stop("SampleSpecies must be one of the following: Human, Mouse.")
}
TopGene <- c()
TopNumber <- 10

Idents(scObject) <- paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")
for (i in levels(scObject)) {
  CalculateResult <- CalculateMarker(object = scObject,
                                     ClusterName = paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony"),
                                     ClusterName.1 = i,
                                     min.pct = 0.01, logfc.threshold = 0.01,
                                     max.cells.per.ident = 3000, only.pos = TRUE,
                                     annotation = AnnoList)
  write.csv(CalculateResult, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Marker.csv"),
            quote = TRUE, row.names = FALSE)
  
  TopGene <- c(TopGene, top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName)
  
  {
    p1 <- VlnPlot(scObject, features = top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName,
                  pt.size = 0.01, assay = "RNA", raster = FALSE) &
      NoLegend() & labs(x = NULL) & scale_fill_manual(values = PlotTheme$Color) &
      PlotTheme$Box & PlotTheme$FontColor
    ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_ViolinPlot.png"),
           bg = "white", width = 17, height = 10)
    ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_ViolinPlot.pdf"),
           bg = "white", width = 17, height = 10)
    
    p1 <- FeaturePlot(scObject, features = top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName, reduction = "umap",
                      raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),
                      pt.size = 0.05, max.cutoff = "q99") &
      PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
    ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_FeaturePlot.png"),
           bg = "white", width = 16, height = 10)
    ggsave(p1, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Top", TopNumber, "_FeaturePlot.pdf"),
           bg = "white", width = 16, height = 10)
  }
  
  CalculateResult <- filter(CalculateResult, p_val < 0.05)
  CalculateResult <- filter(CalculateResult, avg_log2FC >= 0.26)
  write.csv(CalculateResult, file = paste0(ResultPath, "res.", sprintf("%0.1f", ResNubmer), "_Cluster", i, "_Marker_pvalue_0.05_log2FC_0.26.csv"),
            quote = TRUE, row.names = FALSE)
}
## CelltypeMarker -----
suppressMessages({
  library(dplyr)
  library(stringr)
})

source("/home/kaicheng/fuction/CalculateMarker.R")

ResultPath <- c("./CelltypeMarker/1.Marker/")

##
{
  ## Create result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  
  ##
  Idents(scObject) <- "Celltype"
  if (str_to_title(SampleSpecies) == "Human") {
    AnnoList <- read.table("/home/kaicheng/file/Human_Ensembl_Symbol_ENTREZID.xls", header = TRUE, sep = "\t", quote = "")
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    AnnoList <- read.table("/home/kaicheng/file/Mouse_Ensembl_Symbol_ENTREZID.xls", header = TRUE, sep = "\t", quote = "")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  TopGene <- c()
  TopNumber <- 10
  
  for (i in levels(scObject)) {
    CalculateResult <- CalculateMarker(object = scObject,
                                       ClusterName = "Celltype", ClusterName.1 = i,
                                       min.pct = 0.01, logfc.threshold = 0.01,
                                       max.cells.per.ident = 3000, only.pos = TRUE,
                                       annotation = AnnoList)
    write.csv(CalculateResult, file = paste0(ResultPath, i, "_Marker.csv"),
              quote = TRUE, row.names = FALSE)
    
    TopGene <- c(TopGene, top_n(CalculateResult, n = TopNumber, wt = avg_log2FC)$GeneName)
    
    CalculateResult <- filter(CalculateResult, p_val < 0.05)
    CalculateResult <- filter(CalculateResult, avg_log2FC >= 0.26)
    write.csv(CalculateResult, file = paste0(ResultPath, i, "_Marker_pvalue_0.05_log2FC_0.26.csv"),
              quote = TRUE, row.names = FALSE)
  }
}

{
  {
    p1 <- DotPlot(scObject, features = unique(TopGene), cols = c("white", "red"),
                  cluster.idents = FALSE, assay = "RNA", dot.scale = 4, col.min = 0) &
      labs(x = NULL, y = NULL) & PlotTheme$Box & PlotTheme$FontColor &
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    p1[["data"]][["id"]] <- factor(p1[["data"]][["id"]], levels = rev(levels(p1[["data"]][["id"]])))
    ggsave(p1, file = paste0(ResultPath, "CelltypeMarker_Top", TopNumber, "_Dotplot.png"),
           bg = "white", width = 15, height = 5)
    ggsave(p1, file = paste0(ResultPath, "CelltypeMarker_Top", TopNumber, "_Dotplot.pdf"),
           bg = "white", width = 15, height = 5)
  }
  
  {
    cluster.averages <- as.data.frame(t(AverageExpression(scObject)[["RNA"]][unique(TopGene),]))
    
    p1 <- pheatmap::pheatmap(cluster.averages, scale = "column", cluster_rows = FALSE, cluster_cols = FALSE,
                             border_color = NA, angle_col = 90)
    ggsave(p1, file = paste0(ResultPath, "CelltypeMarker_Top", TopNumber, "_Heatmap.png"),
           bg = "white", width = 15, height = 5)
    ggsave(p1, file = paste0(ResultPath, "CelltypeMarker_Top", TopNumber, "_Heatmap.pdf"),
           bg = "white", width = 15, height = 5)
  }
  
  {
    NeedGene <- list(B = c("CD79A", "MS4A1", "IGHD"),
                     `T` = c("CD3D", "CD3E", "CD8A"),
                     Myeloid = c("LYZ", "CD14", "FCGR3A"),
                     Neutrophil = c("CSF3R", "S100A8", "S100A9"),
                     Mast = c("KIT", "CPA3", "TPSAB1"),
                     Medulla = c("EPCAM", "CHGA", "TH"),
                     Fibroblast = c("COL1A1", "ACTA2", "DCN"),
                     Endothelial = c("PECAM1", "CLDN5", "VWF"),
                     Neuron = c("NRXN1", "NRXN3", "CDH19"))
    p1 <- DotPlot(scObject, features = NeedGene, cols = c("white", "red"),
                  cluster.idents = FALSE, assay = "RNA", dot.scale = 7, col.min = 0) &
      labs(x = NULL, y = NULL) & PlotTheme$Box & PlotTheme$FontColor &
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
    p1[["data"]][["id"]] <- factor(p1[["data"]][["id"]], levels = rev(levels(p1[["data"]][["id"]])))
    ggsave(p1, file = paste0(ResultPath, "Celltype_CanonicalMarker_Dotplot.png"),
           bg = "white", width = 20, height = 6)
    ggsave(p1, file = paste0(ResultPath, "Celltype_CanonicalMarker_Dotplot.pdf"),
           bg = "white", width = 20, height = 6)
  }
}

## Celltype_GO_KEGG -----
suppressMessages({
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
})

source("../../Function/RunGO.R")
source("../../Function/RunKEGG.R")

CalculatePath <- c("./CelltypeMarker/1.Marker/")
ResultPathGO <- c("./CelltypeMarker/2.GO/")
ResultPathKEGG <- c("./CelltypeMarker/3.KEGG/")

##
{
  ## Get the file name containing only significant information
  FileName <- list.files(CalculatePath)
  FileName <- FileName[str_detect(FileName, ".csv")]
  FileName <- FileName[str_detect(FileName, "pvalue", negate = FALSE)]
  FileName <- gsub(".csv", "", FileName)
  
  ## Create GO/KEGG result path
  dir.create(ResultPathGO, showWarnings = FALSE, recursive = TRUE)
  dir.create(ResultPathKEGG, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./Rds/", ResultPathGO), showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./Rds/", ResultPathKEGG), showWarnings = FALSE, recursive = TRUE)
  
  if (str_to_title(SampleSpecies) == "Human") {
    suppressMessages({
      library(org.Hs.eg.db)
    })
    
    NeedOrgDb <- org.Hs.eg.db
    NeedOrganism <- "hsa"
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    suppressMessages({
      library(org.Mm.eg.db)
    })
    
    NeedOrgDb <- org.Mm.eg.db
    NeedOrganism <- "mmu"
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  ## Run GO
  RunGO(FileName = FileName,
        PDF = OUTPUTPDF,
        FilePath = CalculatePath,
        ResultPath = ResultPathGO,
        OrgDb = NeedOrgDb,
        UpGene = TRUE, DownGene = FALSE)
  
  ## Run KEGG
  RunKEGG(FileName = FileName,
          PDF = OUTPUTPDF,
          FilePath = CalculatePath,
          ResultPath = ResultPathKEGG,
          organism = NeedOrganism,
          UpGene = TRUE, DownGene = FALSE)
}

## Celltype_GSEA -----
suppressMessages({
  library(clusterProfiler)
  library(dplyr)
  library(stringr)
})

source("../../Function/RunGSEA.R")

CalculatePath <- c("./CelltypeMarker/1.Marker/")
ResultPath <- c("./CelltypeMarker/4.GSEA/")

##
{
  ## Get the file name containing all information
  FileName <- list.files(CalculatePath)
  FileName <- FileName[str_detect(FileName, ".csv")]
  FileName <- FileName[str_detect(FileName, "pvalue", negate = TRUE)]
  FileName <- gsub(".csv", "", FileName)
  
  ## Create GESA result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)
  dir.create(paste0("./Rds/", ResultPath), showWarnings = FALSE, recursive = TRUE)
  
  if (str_to_title(SampleSpecies) == "Human" | str_to_title(SampleSpecies) == "Mouse") {
    ## Run GESA
    RunGSEA(FileName = FileName,
            PDF = OUTPUTPDF,
            FilePath = CalculatePath,
            ResultPath = ResultPath,
            organism = str_to_title(SampleSpecies),
            GESA_H = TRUE, GESA_C2 = TRUE, GESA_C3 = TRUE,
            GESA_C5 = TRUE, GESA_C8 = TRUE)
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
}