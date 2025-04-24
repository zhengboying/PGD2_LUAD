suppressMessages({
  library(dplyr)
  library(stringr)
  library(reshape2)
  library(Seurat)
  library(ggplot2)
})

## setup
{
  
    project_path =  "/Users/zhengboying/Documents/DMPJ0001-PGD2-byzheng/Analysis/GSE131907/"
    setwd(project_path)
    source("../../Function/PlotTheme.R")
  
  ResultPath <- c("./4.Subset/") 
}
##############################subset obj ----------------
if(TRUE){
  scObject<-readRDS("./Rds/Seurat/4.umap_tsne_nFeature2000.Rds")
  #24636 features across 185197 samples within 1 assay 
  Idents(scObject) <- "RNA_snn_res.0.1Harmony"
  
}


#remove cluster5 4337 cell mix mutiple cluster 
cluster5_cells <- subset(scObject, idents = 5)
write.csv(cluster5_cells@meta.data, file = paste0(ResultPath,"filteredcell/cluster5_cells.csv"))

scObject<-subset(scObject, idents = c(0,1,2,3,4,6,7,8,9,10,11,12,13,14))
##24636 features across 180860 samples within 1 assay


## 剔除cluster3的CD3D+细胞
{
  cluster3_cells <- subset(scObject, idents = 3)
  cluster3_cells_mixT <- subset(cluster3_cells, subset = (CD3D > 1.5) | (CD3E > 1.5)|(PTPRC > 1.5)) ##3818
  write.csv(cluster3_cells_mixT@meta.data, file = paste0(ResultPath,"filteredcell/cluster3_cells.csv"))
  scObject <- subset(scObject, cells = setdiff(colnames(scObject), colnames(cluster3_cells_mixT)))
  #24636 features across 177042 samples within 1 assay
}

{## check
  reduction = "UMAP"
    CanonicalMarker<-c("EPCAM","KRT19","CDH1",
                       "PTPRC","MKI67",
                       "LYZ","CD68",
                       "CD3D","CD3E")
  
  for (i in CanonicalMarker) {
    try({
      p1 <- FeaturePlot(scObject, features = i, reduction = reduction,
                        raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),
                        pt.size = 0.5, max.cutoff = "q99") &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot.png"),
             bg = "white", width = 8, height = 7)
      if (plot_density == TRUE){
        p1 <- plot_density(scObject, features = i, reduction = reduction, size = 1) &
          PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
        ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot_Density.png"),
               bg = "white", width = 8, height = 7)
      }
    })
  }
}

{
  ##Recreate SeuratObject
  scObject <- CreateSeuratObject(counts = scObject@assays$RNA$counts,
                                   min.cells = 3, min.features = 1,
                                   meta.data = FetchData(scObject, vars = c("Sample")))
    Idents(scObject) <- "Sample"
    scObject$Sample <- scObject@active.ident
  }
  saveRDS(scObject, file = paste0("./Rds/Seurat/2.1_filtercells.Rds"))



# 06.Basic analysis -----



suppressMessages({
  library(stringr)
})

## Setup
{
  ## FindVariableFeatures Parameters
  SampleSpecies = "Human"
  Test_nFeature =FALSE
  if (Test_nFeature == TRUE){
    FindVariableFeatures_nfeatures <- c(1500,2000,3000,5000)
  }else{
    FindVariableFeatures_nfeatures <- c(3000)
  }
  
  FindVariableFeatures_selection.method = "vst"
  
  ## ScaleData Parameters
  ScaleData_vars.to.regress = c("nCount_RNA", "nFeature_RNA", "CC.Difference")
  
  ## PCA Parameters
  PCA_npcs = 30
}

## 06.1 CellCycle (False)-----
if(TRUE){
  
  if (str_to_title(SampleSpecies) == "Human") {
    load("../../Database/Human_CellCycleGene.rda")
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    load("../../Database/Mouse_CellCycleGene.rda")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  scObject <- NormalizeData(scObject, assay = "RNA", verbose = FALSE)
  scObject <- CellCycleScoring(scObject, set.ident = TRUE,
                               s.features = rownames(scObject)[rownames(scObject) %in% CellCycleGene$s.genes],
                               g2m.features = rownames(scObject)[rownames(scObject) %in% CellCycleGene$g2m.genes])
  
  ## CC.Difference is the difference between the S and G2M scores
  scObject$CC.Difference <- scObject$S.Score - scObject$G2M.Score
  
  Idents(scObject) <- "Phase"
  Idents(scObject) <- factor(Idents(scObject), levels = sort(levels(scObject)))
  scObject$Phase <- scObject@active.ident
}

## 06.2 Basic Seurat -----
for (n in FindVariableFeatures_nfeatures){
  ## FindVariableFeatures
  scObject <- FindVariableFeatures(scObject, assay = "RNA", selection.method = FindVariableFeatures_selection.method,
                                   nfeatures = FindVariableFeatures_nfeatures, verbose = FALSE)
  
  ## Remove bad genes
  if (str_to_title(SampleSpecies) == "Human") {
    VariableFeatures(scObject) <- grep(pattern = "^MT-", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^RP[SL]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^HB[ABD]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^IG[HKL]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "AC233755.1", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^TR[ABDG]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    VariableFeatures(scObject) <- grep(pattern = "^mt-", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Rp[sl]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Hb[abd]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Ig[hkl]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
    VariableFeatures(scObject) <- grep(pattern = "^Tr[abdg]", VariableFeatures(scObject), value = TRUE, invert = TRUE)
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
  
  ## ScaleData
  scObject <- ScaleData(scObject, assay = "RNA", features = VariableFeatures(scObject),
                        vars.to.regress = ScaleData_vars.to.regress,
                        verbose = FALSE)
  
  ## PCA
  scObject <- RunPCA(scObject, assay = "RNA", npcs = PCA_npcs,
                     seed.use = 16, verbose = FALSE)
  
  ## Save PCA result
  saveRDS(scObject, file = paste0("./Rds/Seurat/3.1_filtercells_PCA_nFeature",n,".Rds"))
}



# 07. Cluster -----
if (TRUE){
  scObject<-readRDS("./Rds/Seurat/3.1_filtercells_PCA_nFeature3000.Rds")
}
suppressMessages({
  library(harmony)
  library(dplyr)
  library(reshape2)
})

## Setup
{
  source("../../Function/ClusterStat.R")
  source("../../Function/ClusterStatPlot.R")
  source("../../Function/PlotTheme.R")
  ResultPath <- c("./5.reCluster/")
  HarmonyIntegration = TRUE
  NoIntegration = FALSE
  OUTPUTPDF = FALSE
  ## NoIntegration and HarmonyIntegration parameters
  dims = 1:30
  RNA_snn_res = c( 0.1,0.4,0.6)
  TSNE = FALSE
  UMAP = TRUE
  min.dist = 0.5
  NoIntegration = FALSE
  ## Load data
  nFeature = 3000
  #scObject<-readRDS(paste0("./Rds/Seurat/3.PCA_nFeature",nFeature,".Rds"))
}

## 07.1 NoIntegration -----
if (NoIntegration == TRUE){
  ## Create result path
  dir.create(paste0(ResultPath,"NoIntegration/"), showWarnings = FALSE, recursive = TRUE)
  
  ## RUN UMAP and tSNE
  scObject <- FindNeighbors(scObject, reduction = "pca", dims = dims,
                            verbose = FALSE)
  if(UMAP == TRUE){
    scObject <- RunUMAP(scObject, assay = "RNA", reduction = "pca", dims = dims,
                        n.components = 2, n.neighbors = 40, min.dist = min.dist,
                        reduction.name = "umap", reduction.key = "umap_",
                        seed.use = 16, verbose = FALSE)
  }
  
  if(TSNE == TRUE){
    scObject <- RunTSNE(scObject, assay = "RNA", reduction = "pca", dims = dims,
                        dim.embed = 2,
                        reduction.name = "tSNE", reduction.key = "tSNE_",
                        seed.use = 16, verbose = FALSE)
  }
  
  
  ## Plot mix umap and tsne
  for (i in c("Sample", "Phase")) {
    if(UMAP == TRUE){
      p1 <- DimPlot(scObject, reduction = "umap", group.by = i,
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
        theme(legend.title = element_blank())
      if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
      ggsave(p1, file = paste0(ResultPath, "NoIntegration/",i, "_umap.png"),
             bg = "white", width = 8, height = 7)
    }
    
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = "tSNE", group.by = i,
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
        theme(legend.title = element_blank())
      if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
      ggsave(p1, file = paste0(ResultPath, "NoIntegration/",i, "_tsne.png"),
             bg = "white", width = 8, height = 7)
    }
  }
  
  ##
  for(i in RNA_snn_res) {
    scObject <- FindClusters(scObject, resolution = i, random.seed = 16, verbose = FALSE)
    
    scObject@meta.data[[paste0("RNA_snn_res.", i)]] <- NULL
    scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW")]] <- scObject$seurat_clusters
    scObject$seurat_clusters <- NULL
    
    ## Plot
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = "umap", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"),
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "NoIntegration/","res.", sprintf("%0.1f", i), "_umap.png"),
             bg = "white", width = 8, height = 7)
      
    }
    if(TSNE == TRUE){
      p1 <- DimPlot(scObject, reduction = "tSNE", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"),
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath,"NoIntegration/", "res.", sprintf("%0.1f", i), "_tsne.png"),
             bg = "white", width = 8, height = 7)
    }
    
    ## Calculate number
    ClusterStatResult <- ClusterStat(object = scObject, percent = FALSE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "NoIntegration/", "res.", sprintf("%0.1f", i), "_Sample_Stat.csv"),
              quote = TRUE, row.names = FALSE)
    
    ## Plot base on number
    ClusterStatPlot(PlotData = ClusterStatResult,PDF = OUTPUTPDF,
                    RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"), ColType = "Sample",
                    PlotPrefixName = paste0("res.", sprintf("%0.1f", i)),
                    RowTypeAngle = FALSE, ColTypeAngle = TRUE)
    
    ## Calculate percent
    ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "RAW"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath,"NoIntegration/", "res.", sprintf("%0.1f", i), "_Sample_Stat_Percent.csv"),
              quote = TRUE, row.names = FALSE)
  }
}

## 07.2 HarmonyIntegration -----
if (HarmonyIntegration == TRUE){
  
  ## Create result path
  dir.create(paste0(ResultPath,"HarmonyIntegration/"), showWarnings = FALSE, recursive = TRUE)
  
  ##
  scObject <- RunHarmony(scObject, assay.use = "RNA", group.by.vars = "Sample",
                         theta = 1, lambda = 0.5, max.iter = 20,
                         early_stop = -Inf, verbose = FALSE) # epsilon.harmony>early_stop
  
  scObject <- FindNeighbors(scObject, reduction = "harmony", dims = dims,
                            verbose = FALSE)
  if(UMAP == TRUE){
    scObject <- RunUMAP(scObject, assay = "RNA", reduction = "harmony", dims = dims,
                        n.components = 2, n.neighbors = 40, min.dist = min.dist,
                        reduction.name = "UMAP", reduction.key = "UMAP_",
                        seed.use = 16, verbose = FALSE)
  }
  
  if(TSNE == TRUE){
    scObject <- RunTSNE(scObject, assay = "RNA", reduction = "harmony", dims = dims,
                        dim.embed = 2,
                        reduction.name = "TSNE", reduction.key = "TSNE_",
                        seed.use = 16, verbose = FALSE)
  }
  
  
  ## Plot mix umap and tsne
  
  for (i in c("Sample", "Phase")) {
    if(UMAP == TRUE){
      p1 <- DimPlot(scObject, reduction = "UMAP", group.by = i,
                    label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
        theme(legend.title = element_blank())
      if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
      ggsave(p1, file = paste0(ResultPath,"HarmonyIntegration/", i, "_umap.png"),
             bg = "white", width = 8, height = 7)
    }
    
    if (TSNE == TRUE) {
      p1 <- DimPlot(scObject, reduction = "TSNE", group.by = i,
                    label = FALSE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        labs(title = NULL) & PlotTheme$Box & PlotTheme$FontColor &
        theme(legend.title = element_blank())
      if (i != "Phase") { p1 <- p1 & scale_color_manual(values = PlotTheme$Color) }
      ggsave(p1, file = paste0(ResultPath,"HarmonyIntegration/", i, "_tsne.png"),
             bg = "white", width = 8, height = 7)
    }
    
  }
  
  ## Stat and plot
  for(i in RNA_snn_res) {
    scObject <- FindClusters(scObject, resolution = i, random.seed = 16, 
                             algorithm = 1, n.start = 10, n.iter = 10,verbose = FALSE)
    
    scObject@meta.data[[paste0("RNA_snn_res.", i)]] <- NULL
    scObject@meta.data[[paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony")]] <- scObject$seurat_clusters
    scObject$seurat_clusters <- NULL
    
    ## Plot
    if(UMAP == TRUE){
      p1 <- DimPlot(scObject, reduction = "UMAP", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"),
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath,"HarmonyIntegration/", "res.", sprintf("%0.1f", i), "_umap.png"),
             bg = "white", width = 8, height = 7)
    }
    
    if (TSNE == TRUE) {
      p1 <- DimPlot(scObject, reduction = "TSNE", group.by = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"),
                    label = TRUE, shuffle = TRUE, raster = FALSE, repel = FALSE,
                    pt.size = 0.5, label.size = 5, seed = 16) &
        scale_color_manual(values = PlotTheme$Color) &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "HarmonyIntegration/","res.", sprintf("%0.1f", i), "_tsne.png"),
             bg = "white", width = 8, height = 7)
    }
    
    
    ## Calculate number
    ClusterStatResult <- ClusterStat(object = scObject, percent = FALSE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "HarmonyIntegration/","res.", sprintf("%0.1f", i), "_Sample_Stat.csv"),
              quote = TRUE, row.names = FALSE)
    
    ## Plot base on number
    ClusterStatPlot(PlotData = ClusterStatResult,PDF = OUTPUTPDF,
                    RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"), ColType = "Sample",
                    PlotPrefixName = paste0("res.", sprintf("%0.1f", i)),
                    RowTypeAngle = FALSE, ColTypeAngle = TRUE,
                    ResultPath = paste0(ResultPath, "HarmonyIntegration/"))
    
    
    ## Calculate percent
    ClusterStatResult <- ClusterStat(object = scObject, percent = TRUE,
                                     RowType = paste0("RNA_snn_res.", sprintf("%0.1f", i), "Harmony"), ColType = "Sample")
    write.csv(ClusterStatResult, file = paste0(ResultPath, "HarmonyIntegration/","res.", sprintf("%0.1f", i), "_Sample_Stat_Percent.csv"),
              quote = TRUE, row.names = FALSE)
  }
  
}

## Save step 07 result
{
  saveRDS(scObject, file = paste0("./Rds/Seurat/4.1_filtercells_nFeature",nFeature,".Rds"))
}



# 09.CanonicalMarker -----

suppressMessages({
  library(Nebulosa)##plot_density
  library(stringr)
  library(RColorBrewer)
  library(tidyverse)
  library(ComplexHeatmap)
})

## Setup
{
  source("../../Function/Plot.R")
  ResultPath <- c("./3.res.Marker/")
  FeaturePlot = TRUE
  DotPlot = TRUE ## surpport Human only.
  plot_density = FALSE ## plot_density耗时较久，非交付需求可选择FALSE
  ResNubmer = 0.1 ## from step 08
  reduction = "UMAP" ##tSNE,UMAP,umap
}  

## Load CanonicalMarker  
{
  ## Create result path
  dir.create(ResultPath, showWarnings = FALSE, recursive = TRUE)  
  ##
  if (str_to_title(SampleSpecies) == "Human") {
    ## all markers,which need to be filtered after check results,
    CanonicalMarker<-read.csv("../../Database/CellType.annotdata/ref_cell_type.csv",header = T)
  } else if (str_to_title(SampleSpecies) == "Mouse") {
    source("../../Database/Mouse_CanonicalMarker.R")
  } else {
    stop("SampleSpecies must be one of the following: Human, Mouse.")
  }
}  

## 09.1 CanonicalMarker DotPlot -----
{
  if (str_to_title(SampleSpecies) == "Human") {
    col_plates_init<-c("NK" = "grey","T"="#90E663","Leukocyte" = "#39A093",
                       "B"="#90E663","Plasma"="#11E693","Oligo" = "brown",
                       "Mast" = "#63B3E6","Myeloid" = "#13F4FF","Monocyte" = "#78C4FF", 
                       "pDC" = "#13F4FF","DC" = "#78C4FF","Neutrophil" = "pink","Pericyte"= "white",
                       "Epithelium" = "#AB49F0", "Platelet" = "black","Megakaryocyte" = "orange",
                       "Endothelial"="#F5802C","MSCs" = "#E5202C","Stromal" = "#F5E02C","Fibroblast" = "#DEF02C")
    
    Idents(scObject)<-paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony")
    ht<-anno_hematmap(obj = scObject,
                      markers = CanonicalMarker,
                      col_plates = col_plates_init,
                      cell_cluster = 9,#not used
                      group = paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony"),
                      row_names_fontsize = 8,
                      seed = 12)
    ##save
    png(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot.png"), width =12, height = 16, units = "in", res = 600)
    draw(ht)
    dev.off()  
    pdf(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot.pdf"), width =12, height = 16)
    draw(ht)
    dev.off()
    
    ## Filtered markers
    if(TRUE){
      #去掉CellTye =neutrophils,Pericytes,mregDCs,Lam,"RTMs",CMono"
      #CanonicalMarker<-CanonicalMarker[!CanonicalMarker$CellType %in% c("Pericytes","DC"),]#
      ##删掉一些不怎么表达的基因
      CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CMA1","KRT1","IL1RAPL1","PF4","FOLR2","FABP3",
                                                                      "PLA2G2D","SH2D1B","RGS5","IGLL1","AMPD1","PPBP","TBX21",
                                                                      "TNFRSF4","BATF","FOXP3","CTLA4","LAIR2","EPCAP"),]
      
      CanonicalMarker<-CanonicalMarker[!CanonicalMarker$Marker %in% c("CD34","CDH5","CLEC9A","CLEC10A","CD1C","IL1RL1",
                                                                      "CTSG","HDC","CPA3","S100A8","S100A9","COL1A1","COL3A1"),]
      
      ht<-anno_hematmap(obj = scObject,
                        markers = CanonicalMarker,
                        col_plates = col_plates_init,
                        cell_cluster = 9,#not used
                        group = paste0("RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony"),
                        row_names_fontsize = 8,
                        seed = 12)
      ##save
      png(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot_Filterd.png"), width =12, height = 16, units = "in", res = 600)
      draw(ht)
      dev.off()  
      pdf(file = paste0(ResultPath, "RNA_snn_res.", sprintf("%0.1f", ResNubmer), "Harmony_","CanonicalMarkerDotPlot_Filterd.pdf"), width =12, height = 16)
      draw(ht)
      dev.off()
      
    }
  }
}

## 09.2 CanonicalMarker FeaturePlot ----
{
  ##
  if (str_to_title(SampleSpecies) == "Human") {
  for (i in CanonicalMarker) {
    try({
      p1 <- FeaturePlot(scObject, features = i, reduction = reduction,
                        raster = FALSE, order = TRUE, cols = c("lightgrey", "red"),
                        pt.size = 0.5, max.cutoff = "q99") &
        PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
      ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot.png"),
             bg = "white", width = 8, height = 7)
      if (plot_density == TRUE){
        p1 <- plot_density(scObject, features = i, reduction = reduction, size = 1) &
          PlotTheme$Box & PlotTheme$FontColor & theme(legend.title = element_blank())
        ggsave(p1, file = paste0(ResultPath, "CanonicalMarker/", i, "_FeaturePlot_Density.png"),
               bg = "white", width = 8, height = 7)
      }
    })
  }
}
